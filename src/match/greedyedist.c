/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <limits.h>
#include <stdio.h>
#include "core/chardef.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/types_api.h"
#include "greedyedist.h"

#define GT_FRONT_STORE(GL,F,V)         F = V
#define GT_FRONT_ROWVALUE(FVAL)        *(FVAL)
#define GT_FRONT_MINUSINFINITY(GL)     ((GL)->integermin)

typedef long GtFrontvalue;

typedef struct
{
  long offset,   /* absolute base address of the front */
       width,    /* width of the front */
       left;     /* left boundary (negative), */
                 /* -left is relative address of entry */
} GtFrontspec;

/*
  A structure to store the global values.
*/

struct GtFrontResource
{
  unsigned long currentallocated;
  long ulen,
       vlen,
       integermin;
  GtFrontvalue *frontspace;
};

static unsigned long sumofoddnumbers(unsigned long max)
{
  return max * max;
}

GtFrontResource *gt_frontresource_new(unsigned long maxdist)
{
  GtFrontResource *ftres = gt_malloc(sizeof *ftres);

  ftres->currentallocated = sumofoddnumbers(maxdist);
  ftres->frontspace = gt_malloc(sizeof *ftres->frontspace
                                * ftres->currentallocated);
  return ftres;
}

void gt_frontresource_delete(GtFrontResource *ftres)
{
  if (ftres != NULL)
  {
    gt_free(ftres->frontspace);
    gt_free(ftres);
  }
}

#ifdef SKDEBUG
static void showfront(const GtFrontResource *ftres,
                      const GtFrontspec *fspec,
                      long r)
{
  GtFrontvalue *fval;
  long k;

  for (fval = ftres->frontspace + fspec->offset, k=fspec->left;
       k < fspec->left + fspec->width; k++, fval++)
  {
    if (r <= 0 || k <= -r || k >= r)
    {
      printf("(k=%ld)%ld ",k,GT_FRONT_ROWVALUE(fval));
    } else
    {
      printf("(k=%ld)undef ",k);
    }
  }
  (void) putchar('\n');
}
#endif

static void frontspecparms(const GtFrontResource *ftres,
                           GtFrontspec *fspec,
                           long p,
                           long r)
{
  if (r <= 0)
  {
    fspec->left = -p;
    fspec->width = p + p + 1;
  } else
  {
    fspec->left = MAX(-ftres->ulen,-p);
    fspec->width = MIN(ftres->vlen,p) - fspec->left + 1;
  }
#ifdef SKDEBUG
  printf("p=%ld,offset=%ld,left=%ld,width=%ld\n",p,
                                                 fspec->offset,
                                                 fspec->left,
                                                 fspec->width);
#endif
}

static long accessfront(const GtFrontResource *ftres,
                        GtFrontvalue *fval,
                        const GtFrontspec *fspec,
                        long k)
{
  if (fspec->left <= k && k < fspec->left + fspec->width)
  {
    return GT_FRONT_ROWVALUE(fval+k);
  }
  return GT_FRONT_MINUSINFINITY(ftres);
}

/*
  The following function evaluates an entry \(front(k,p)\) in a
  forward direction.
*/

static void evalentryforward(const GtSeqabstract *useq,
                             const GtSeqabstract *vseq,
                             GtFrontResource *ftres,
                             GtFrontvalue *fval,
                             const GtFrontspec *fspec,
                             long k)
{
  long value, t;
  GtFrontvalue *fptr;

#ifdef SKDEBUG
  printf("evalentryforward(k=%ld)\n",k);
#endif
  fptr = ftres->frontspace + fspec->offset - fspec->left;
  t = accessfront(ftres,fptr,fspec,k) + 1;         /* same diagonal */
#ifdef SKDEBUG
  printf("same: access(k=%ld)=%ld\n",k,t-1);
#endif

  value = accessfront(ftres,fptr,fspec,k-1);       /* diagonal below */
#ifdef SKDEBUG
  printf("below: access(k=%ld)=%ld\n",k-1,value);
#endif
  if (t < value)
  {
    t = value;
  }
  value = accessfront(ftres,fptr,fspec,k+1) + 1;     /* diagonal above */
#ifdef SKDEBUG
  printf("above: access(k=%ld)=%ld\n",k+1,value-1);
#endif
  if (t < value)
  {
    t = value;
  }
#ifdef SKDEBUG
  printf("maximum: t=%ld\n",t);   /* the maximum over three values */
#endif
  if (t < 0 || t+k < 0)             /* no negative value */
  {
    GT_FRONT_STORE(ftres,GT_FRONT_ROWVALUE(fval),GT_FRONT_MINUSINFINITY(ftres));
  } else
  {
    if (ftres->ulen != 0 && ftres->vlen != 0 &&
        t < ftres->ulen && t + k < ftres->vlen)
    {
      unsigned long lcp, minlen
        = (unsigned long) MIN(ftres->ulen - t,ftres->vlen - (t + k));
      lcp = gt_seqabstract_lcp(true,
                               useq,
                               vseq,
                               (unsigned long) t,
                               (unsigned long) (t + k),
                               minlen);
      t += lcp;
    }
    if (t > ftres->ulen || t + k > ftres->vlen)
    {
      GT_FRONT_STORE(ftres,GT_FRONT_ROWVALUE(fval),
                     GT_FRONT_MINUSINFINITY(ftres));
    } else
    {
      GT_FRONT_STORE(ftres,GT_FRONT_ROWVALUE(fval),t);
    }
  }
}

/*
  The following function evaluates a front in forward direction.
  It returns true if any of the returned values is at least 0.
*/

static bool evalfrontforward(const GtSeqabstract *useq,
                             const GtSeqabstract *vseq,
                             GtFrontResource *ftres,
                             const GtFrontspec *prevfspec,
                             const GtFrontspec *fspec,
                             long r)
{
  long k;
  bool defined = false;
  GtFrontvalue *fval;

  for (fval = ftres->frontspace + fspec->offset, k = fspec->left;
       k < fspec->left + fspec->width; k++, fval++)
  {
    if (r <= 0 || k <= -r || k >= r)
    {
      evalentryforward(useq,vseq,ftres,fval,prevfspec,k);
      if (GT_FRONT_ROWVALUE(fval) >= 0)
      {
        defined = true;
      }
#ifdef SKDEBUG
      printf("store front[k=%ld]=%ld ",k,GT_FRONT_ROWVALUE(fval));
      printf("at index %ld\n",(long) (fval-ftres->frontspace));
#endif
    } else
    {
#ifdef SKDEBUG
      printf("store front[k=%ld]=GT_FRONT_MINUSINFINITY ",k);
      printf("at index %ld\n",(long) (fval-ftres->frontspace));
#endif
      GT_FRONT_STORE(ftres,GT_FRONT_ROWVALUE(fval),
                     GT_FRONT_MINUSINFINITY(ftres));
    }
  }
#ifdef SKDEBUG
  printf("frontvalues[r=%ld]=",r);
  showfront(ftres,fspec,r);
#endif
  return defined;
}

/*
  The following function evaluates the entry \(front(0,0)\) in a
  forward direction.
*/

static void firstfrontforward(const GtSeqabstract *useq,
                              const GtSeqabstract *vseq,
                              GtFrontResource *ftres,
                              GtFrontspec *fspec)
{
  fspec->left = fspec->offset = 0;
  fspec->width = 1L;
  if (ftres->ulen == 0 || ftres->vlen == 0)
  {
    GT_FRONT_STORE(ftres,GT_FRONT_ROWVALUE(&ftres->frontspace[0]),0);
  } else
  {
    unsigned long lcp = gt_seqabstract_lcp(true,
                                           useq,
                                           vseq,
                                           0,
                                           0,
                                           (unsigned long)
                                           MIN(ftres->ulen,ftres->vlen));
    GT_FRONT_STORE(ftres,GT_FRONT_ROWVALUE(&ftres->frontspace[0]),(long) lcp);
  }
#ifdef SKDEBUG
  printf("forward front[0]=%ld\n",GT_FRONT_ROWVALUE(&ftres->frontspace[0]));
#endif
}

unsigned long greedyunitedist(GtFrontResource *ftres,
                              const GtSeqabstract *useq,
                              const GtSeqabstract *vseq)
{
  unsigned long realdistance, kval;
  long r;
  GtFrontspec frontspecspace[2], *fspec, *prevfspec;
  GtFrontvalue *fptr;

#ifdef SKDEBUG
  printf("unitedistcheckSEPgeneric(ulen=%lu,vlen=%lu)\n",ulenvalue,vlenvalue);
#endif
  gt_assert(gt_seqabstract_length_get(useq) < (unsigned long) LONG_MAX);
  gt_assert(gt_seqabstract_length_get(vseq) < (unsigned long) LONG_MAX);
  ftres->ulen = (long) gt_seqabstract_length_get(useq);
  ftres->vlen = (long) gt_seqabstract_length_get(vseq);
  ftres->integermin = -MAX(ftres->ulen,ftres->vlen);
  prevfspec = &frontspecspace[0];
  firstfrontforward(useq,vseq,ftres,prevfspec);
  if (ftres->ulen == ftres->vlen &&
      GT_FRONT_ROWVALUE(ftres->frontspace) == ftres->vlen)
  {
    realdistance = 0;
  } else
  {
    for (kval=1UL, r=1-MIN(ftres->ulen,ftres->vlen);
         /* Nothing */ ; kval++, r++)
    {
      if (prevfspec == &frontspecspace[0])
      {
        fspec = &frontspecspace[1];
      } else
      {
        fspec = &frontspecspace[0];
      }
      fspec->offset = prevfspec->offset + prevfspec->width;
      frontspecparms(ftres,fspec,(long) kval,r);
      while ((unsigned long) (fspec->offset + fspec->width)
             >= ftres->currentallocated)
      {
        ftres->currentallocated += kval+1;
        ftres->frontspace
          = gt_realloc(ftres->frontspace,
                       sizeof *ftres->frontspace * ftres->currentallocated);
      }
      (void) evalfrontforward(useq,vseq,ftres,prevfspec,fspec,r);
      fptr = ftres->frontspace + fspec->offset - fspec->left;
      if (accessfront(ftres,fptr,fspec,ftres->vlen - ftres->ulen)
          == ftres->ulen)
      {
        realdistance = kval;
        break;
      }
      if (prevfspec == &frontspecspace[0])
      {
        prevfspec = &frontspecspace[1];
      } else
      {
        prevfspec = &frontspecspace[0];
      }
    }
  }
#ifdef SKDEBUG
  printf("unitedistfrontSEP returns %ld\n",realdistance);
#endif
  return realdistance;
}
