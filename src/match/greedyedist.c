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

#define COMPARESYMBOLS(A,B)\
        if ((A) == (GtUchar) SEPARATOR)\
        {\
          gl->ubound = uidx;\
        }\
        if ((B) == (GtUchar) SEPARATOR)\
        {\
          gl->vbound = vidx;\
        }\
        if ((A) != (B) || ISSPECIAL(A))\
        {\
          break;\
        }

#define STOREFRONT(GL,F,V)            F = V
#define ROWVALUE(FVAL)                *(FVAL)

#define MINUSINFINITYFRONT(GL)        ((GL)->integermin)

typedef long Frontvalue;

/*
  A structure to store the global values.
*/

typedef struct
{
  long offset,   /* absolute base address of the front */
       width,    /* width of the front */
       left;     /* left boundary (negative), */
                 /* -left is relative address of entry */
} Frontspec;     /* \Typedef{Frontspec} */

typedef struct
{
  unsigned long ubound,
                vbound;
  long ulen,
       vlen,
       integermin;
  Frontvalue *frontspace;
} FrontResource;

#ifdef SKDEBUG
static void showfront(const FrontResource *gl,
                      const Frontspec *fspec,
                      long r)
{
  Frontvalue *fval;
  long k;

  for (fval = gl->frontspace + fspec->offset, k=fspec->left;
       k < fspec->left + fspec->width; k++, fval++)
  {
    if (r <= 0 || k <= -r || k >= r)
    {
      printf("(k=%ld)%ld ",k,ROWVALUE(fval));
    } else
    {
      printf("(k=%ld)undef ",k);
    }
  }
  (void) putchar('\n');
}
#endif

static void frontspecparms(const FrontResource *gl,
                           Frontspec *fspec,
                           long p,
                           long r)
{
  if (r <= 0)
  {
    fspec->left = -p;
    fspec->width = p + p + 1;
  } else
  {
    fspec->left = MAX(-gl->ulen,-p);
    fspec->width = MIN(gl->vlen,p) - fspec->left + 1;
  }
#ifdef SKDEBUG
  printf("p=%ld,offset=%ld,left=%ld,width=%ld\n",p,
                                                 fspec->offset,
                                                 fspec->left,
                                                 fspec->width);
#endif
}

static long accessfront(const FrontResource *gl,
                        Frontvalue *fval,
                        const Frontspec *fspec,
                        long k)
{
  if (fspec->left <= k && k < fspec->left + fspec->width)
  {
    return ROWVALUE(fval+k);
  }
  return MINUSINFINITYFRONT(gl);
}

/*
  The following function evaluates an entry \(front(k,p)\) in a
  forward direction.
*/

static void evalentryforward(const GtSeqabstract *useq,
                             const GtSeqabstract *vseq,
                             FrontResource *gl,
                             Frontvalue *fval,
                             const Frontspec *fspec,
                             long k)
{
  long value, t;
  unsigned long uidx, vidx;
  GtUchar a, b;
  Frontvalue *fptr;

#ifdef SKDEBUG
  printf("evalentryforward(k=%ld)\n",k);
#endif
  fptr = gl->frontspace + fspec->offset - fspec->left;
  t = accessfront(gl,fptr,fspec,k) + 1;         /* same diagonal */
#ifdef SKDEBUG
  printf("same: access(k=%ld)=%ld\n",k,t-1);
#endif

  value = accessfront(gl,fptr,fspec,k-1);       /* diagonal below */
#ifdef SKDEBUG
  printf("below: access(k=%ld)=%ld\n",k-1,value);
#endif
  if (t < value)
  {
    t = value;
  }
  value = accessfront(gl,fptr,fspec,k+1) + 1;     /* diagonal above */
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
    STOREFRONT(gl,ROWVALUE(fval),MINUSINFINITYFRONT(gl));
  } else
  {
    gt_assert(t >= 0);
    uidx = (unsigned long) t;
    gt_assert(t + k >= 0);
    vidx = (unsigned long) (t + k);
    if (gl->ulen != 0 && gl->vlen != 0)  /* only for nonempty strings */
    {
      for (/* Nothing */; uidx < gl->ubound && vidx < gl->vbound;
           uidx++, vidx++)
      {
        a = gt_seqabstract_encoded_char(useq,uidx);
        b = gt_seqabstract_encoded_char(vseq,vidx);
        COMPARESYMBOLS(a,b);
      }
      t = (long) uidx;
    }
    if (t > (long) gl->ubound || t + k > (long) gl->vbound)
    {
      STOREFRONT(gl,ROWVALUE(fval),MINUSINFINITYFRONT(gl));
    } else
    {
      STOREFRONT(gl,ROWVALUE(fval),t);
    }
  }
}

/*
  The following function evaluates a front in forward direction.
  It returns true if any of the returned values is at least 0.
*/

static bool evalfrontforward(const GtSeqabstract *useq,
                             const GtSeqabstract *vseq,
                             FrontResource *gl,
                             const Frontspec *prevfspec,
                             const Frontspec *fspec,
                             long r)
{
  long k;
  bool defined = false;
  Frontvalue *fval;

  for (fval = gl->frontspace + fspec->offset, k = fspec->left;
       k < fspec->left + fspec->width; k++, fval++)
  {
    if (r <= 0 || k <= -r || k >= r)
    {
      evalentryforward(useq,vseq,gl,fval,prevfspec,k);
      if (ROWVALUE(fval) >= 0)
      {
        defined = true;
      }
#ifdef SKDEBUG
      printf("store front[k=%ld]=%ld ",k,ROWVALUE(fval));
      printf("at index %ld\n",(long) (fval-gl->frontspace));
#endif
    } else
    {
#ifdef SKDEBUG
      printf("store front[k=%ld]=MINUSINFINITYFRONT ",k);
      printf("at index %ld\n",(long) (fval-gl->frontspace));
#endif
      STOREFRONT(gl,ROWVALUE(fval),MINUSINFINITYFRONT(gl));
    }
  }
#ifdef SKDEBUG
  printf("frontvalues[r=%ld]=",r);
  showfront(gl,fspec,r);
#endif
  return defined;
}

/*
  The following function evaluates the entry \(front(0,0)\) in a
  forward direction.
*/

static void firstfrontforward(const GtSeqabstract *useq,
                              const GtSeqabstract *vseq,
                              FrontResource *gl,Frontspec *fspec)
{
  GtUchar a, b;
  unsigned long uidx, vidx;

  fspec->left = fspec->offset = 0;
  fspec->width = 1L;
  if (gl->ulen == 0 || gl->vlen == 0)
  {
    STOREFRONT(gl,ROWVALUE(&gl->frontspace[0]),0);
  } else
  {
    for (uidx = 0, vidx = 0; uidx < gl->ubound && vidx < gl->vbound;
         uidx++, vidx++)
    {
      a = gt_seqabstract_encoded_char(useq,uidx);
      b = gt_seqabstract_encoded_char(vseq,vidx);
      COMPARESYMBOLS(a,b);
    }
    STOREFRONT(gl,ROWVALUE(&gl->frontspace[0]),(long) uidx);
  }
#ifdef SKDEBUG
  printf("forward front[0]=%ld\n",ROWVALUE(&gl->frontspace[0]));
#endif
}

unsigned long greedyunitedist(const GtSeqabstract *useq,
                              const GtSeqabstract *vseq)
{
  unsigned long currentallocated, realdistance;
  FrontResource gl;
  Frontspec frontspecspace[2],
            *fspec,
            *prevfspec;
  Frontvalue *fptr;
  unsigned long kval;
  long r;

#ifdef SKDEBUG
  printf("unitedistcheckSEPgeneric(ulen=%lu,vlen=%lu)\n",ulenvalue,vlenvalue);
#endif
  gt_assert(gt_seqabstract_length_get(useq) < (unsigned long) LONG_MAX);
  gt_assert(gt_seqabstract_length_get(vseq) < (unsigned long) LONG_MAX);
  currentallocated = 1UL;
  gl.frontspace = gt_malloc(sizeof (*gl.frontspace) * currentallocated);
  gl.ubound = gt_seqabstract_length_get(useq);
  gl.vbound = gt_seqabstract_length_get(vseq);
  gl.ulen = (long) gl.ubound;
  gl.vlen = (long) gl.vbound;
  gl.integermin = -MAX(gl.ulen,gl.vlen);
  prevfspec = &frontspecspace[0];
  firstfrontforward(useq,vseq,&gl,prevfspec);
  if (gl.ulen == gl.vlen && ROWVALUE(&gl.frontspace[0]) == gl.vlen)
  {
    realdistance = 0;
  } else
  {
    for (kval=1UL, r=1-MIN(gl.ulen,gl.vlen); /* Nothing */ ; kval++, r++)
    {
      if (prevfspec == &frontspecspace[0])
      {
        fspec = &frontspecspace[1];
      } else
      {
        fspec = &frontspecspace[0];
      }
      fspec->offset = prevfspec->offset + prevfspec->width;
      frontspecparms(&gl,fspec,(long) kval,r);
      while ((unsigned long) (fspec->offset + fspec->width)
             >= currentallocated)
      {
        currentallocated += (kval+1);
        gl.frontspace = gt_realloc(gl.frontspace,sizeof (*gl.frontspace) *
                                   currentallocated);
      }
      (void) evalfrontforward(useq,vseq,&gl,prevfspec,fspec,r);
      fptr = gl.frontspace + fspec->offset - fspec->left;
      if (accessfront(&gl,fptr,fspec,gl.vlen - gl.ulen) == gl.ulen)
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
  gt_free(gl.frontspace);
  return realdistance;
}
