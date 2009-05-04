/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include <stdbool.h>
#include <math.h>
#include "core/ma.h"
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/error_api.h"
#include "core/minmax.h"
#include "divmodmul.h"
#include "intbits-tab.h"
#include "diff-cover.h"
#include "sfx-apfxlen.h"
#include "bcktab.h"
#include "initbasepower.h"
#include "encseq-def.h"

typedef unsigned long Diffvalue;

#define MODV(VAL) ((VAL) & dcov->vmodmask)
#define DIVV(VAL) ((VAL) >> dcov->logmod)

struct Differencecover
{
  unsigned int vparam, logmod, size, vmodmask,
               hvalue,  /* not necessary */
               numofchars,
               prefixlength,
               *coverrank;
  Diffvalue *diffvalues, *diff2pos;
  Seqpos totallength, *sample;
  Bcktab *bcktab;
  const Encodedsequence *encseq;
  unsigned long samplesize;
};

static Diffvalue differencecovertab[] = {
  /*
     1
   */ 0UL,
  /*
     2
   */ 0UL, 1UL,
  /*
     4
   */ 0UL, 1UL, 2UL,
  /*
     8
   */ 0UL, 1UL, 2UL, 4UL,
  /*
     16
   */ 0UL, 1UL, 2UL, 5UL, 8UL,
  /*
     32
   */ 0UL, 1UL, 2UL, 3UL, 7UL, 11UL, 19UL,
  /*
     64
   */ 0UL, 1UL, 2UL, 5UL, 14UL, 16UL, 34UL, 42UL, 59UL,
  /*
     128
   */ 0UL, 1UL, 3UL, 7UL, 17UL, 40UL, 55UL, 64UL, 75UL, 85UL,
    104UL, 109UL, 117UL,
  /*
     256
   */ 0UL, 1UL, 3UL, 7UL, 12UL, 20UL, 30UL, 44UL, 65UL, 80UL,
    89UL, 96UL, 114UL, 122UL, 128UL, 150UL, 196UL, 197UL, 201UL,
    219UL,
  /*
     512
   */ 0UL, 1UL, 2UL, 3UL, 4UL, 9UL, 18UL, 27UL, 36UL, 45UL, 64UL,
    83UL, 102UL, 121UL, 140UL, 159UL, 178UL, 197UL, 216UL, 226UL,
    236UL, 246UL, 256UL, 266UL, 267UL, 268UL, 269UL, 270UL,
  /*
     1024
   */ 0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 6UL, 13UL, 26UL, 39UL, 52UL,
    65UL, 78UL, 91UL, 118UL, 145UL, 172UL, 199UL, 226UL, 253UL,
    280UL, 307UL, 334UL, 361UL, 388UL, 415UL, 442UL, 456UL,
    470UL, 484UL, 498UL, 512UL, 526UL, 540UL, 541UL, 542UL,
    543UL, 544UL, 545UL, 546UL,
  /*
     2048
   */ 0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 6UL, 7UL, 8UL, 9UL, 19UL,
    38UL, 57UL, 76UL, 95UL, 114UL, 133UL, 152UL, 171UL, 190UL,
    229UL, 268UL, 307UL, 346UL, 385UL, 424UL, 463UL, 502UL,
    541UL, 580UL, 619UL, 658UL, 697UL, 736UL, 775UL, 814UL,
    853UL, 892UL, 931UL, 951UL, 971UL, 991UL, 1011UL, 1031UL,
    1051UL, 1071UL, 1091UL, 1111UL, 1131UL, 1132UL, 1133UL,
    1134UL, 1135UL, 1136UL, 1137UL, 1138UL, 1139UL, 1140UL
};

static unsigned int differencecoversizes[]
  = { 1U, 2U, 3U, 4U, 5U, 7U, 9U, 13U, 20U, 28U, 40U, 58U };

static void fillcoverrank(Differencecover *dcov)
{
  unsigned int i, j;

  dcov->coverrank = gt_malloc(sizeof(*dcov->coverrank) * dcov->vparam);
  for (i=0, j=0; i<dcov->vparam; i++)
  {
    dcov->coverrank[i] = j;
    if (j < dcov->size && dcov->diffvalues[j] <= (Diffvalue) i)
    {
      j++;
    }
  }
}

static void filldiff2pos(Differencecover *dcov)
{
  Diffvalue *iptr, *jptr;

  dcov->diff2pos = gt_malloc(sizeof(*dcov->diff2pos) * dcov->vparam);
  for (iptr=dcov->diffvalues + dcov->size - 1; iptr>=dcov->diffvalues; iptr--)
  {
    for (jptr=dcov->diffvalues; jptr<dcov->diffvalues + dcov->size; jptr++)
    {
      dcov->diff2pos[MODV(*jptr - *iptr)] = *iptr;
    }
  }
}

#ifdef WITHcomputehvalue

/* XXX: following function is probably not used */

static unsigned int computehvalue(const Differencecover *dcov,
                                  Seqpos totallength)
{
  Diffvalue next;
  unsigned int h, nmodv = MODV(totallength);

  for (h = 0; h < dcov->size; h++)
  {
    if (dcov->diffvalues[h] <= (Diffvalue) nmodv)
    {
      if (h + 1 < dcov->size)
      {
        next = dcov->diffvalues[h+1];
      } else
      {
        next = (Diffvalue) dcov->vparam;
      }
      if ((Diffvalue) nmodv < next)
      {
        return h;
      }
    }
  }
  gt_assert(false);
  return 0;
}
#endif

Differencecover *differencecover_new(unsigned int vparam,
                                     const Encodedsequence *encseq)
{
  size_t logmod;
  unsigned int offset = 0, v = 1U;
  Differencecover *dcov;
  bool found = false;

  dcov = gt_malloc(sizeof (*dcov));
  for (logmod = 0;
       logmod < sizeof (differencecoversizes)/sizeof (differencecoversizes[0]);
       logmod++)
  {
    if (v == vparam)
    {
      dcov->size = differencecoversizes[logmod];
      dcov->diffvalues = differencecovertab + offset;
      found = true;
      break;
    }
    offset += differencecoversizes[logmod];
    v = MULT2(v);
  }
  if (!found)
  {
    gt_free(dcov);
    return NULL;
  }
  dcov->logmod = (unsigned int) logmod;
  dcov->vparam = 1U << logmod;
  dcov->vmodmask = dcov->vparam-1;
#ifdef WITHcomputehvalue
  dcov->hvalue = computehvalue(dcov,totallength);
#endif
  dcov->sample = NULL;
  dcov->samplesize = 0;
  dcov->totallength = getencseqtotallength(encseq);
  dcov->encseq = encseq;
  dcov->numofchars = getencseqAlphabetnumofchars(encseq);
  fillcoverrank(dcov);
  filldiff2pos(dcov);
  return dcov;
}

unsigned int differencecover_rank(const Differencecover *dcov,Seqpos pos)
{
  return dcov->coverrank[MODV(pos)];
}

unsigned int differencecover_offset(const Differencecover *dcov,
                                    Seqpos pos1,Seqpos pos2)
{
  return (unsigned int) MODV(dcov->diff2pos[MODV(pos2-pos1)] - pos1);
}

void differencecover_delete(Differencecover *dcov)
{
  gt_free(dcov->coverrank);
  gt_free(dcov->diff2pos);
  gt_free(dcov->sample);
  gt_assert(dcov->bcktab != NULL);
  bcktab_delete(&dcov->bcktab);
  gt_free(dcov);
}

unsigned long differencecover_packsamplepos(const Differencecover *dcov,
                                            Seqpos pos)
{
  return dcov->coverrank[MODV(pos)] * (DIVV(dcov->totallength) + 1) +
         (unsigned long) DIVV(pos);
}

static void differencecover_sample(Differencecover *dcov,
                                   Readmode readmode,bool withcheck)
{
  Seqpos pos;
  unsigned int modvalue = 0;
  Diffvalue *diffptr, *afterend;
  unsigned long idx, maxsamplesize;
  Bitsequence *sampleidxused = NULL;
  unsigned int unitsnotspecial;
  Codetype code, *filltable;
  const Codetype **multimappower;
  Encodedsequencescanstate *esr;

  maxsamplesize = (unsigned long) (DIVV(dcov->totallength) + 1) * dcov->size;
  dcov->samplesize = 0;
  dcov->prefixlength = recommendedprefixlength(dcov->numofchars,
                                               (Seqpos) maxsamplesize);
  dcov->bcktab = allocBcktab(dcov->numofchars,
                             dcov->prefixlength,
                             false,
                             NULL,
                             NULL);
  gt_assert(dcov->bcktab != NULL);
  if (withcheck)
  {
    INITBITTAB(sampleidxused,maxsamplesize);
  }
  filltable = filllargestchartable(dcov->numofchars,dcov->prefixlength);
  if (possibletocmpbitwise(dcov->encseq))
  {
    multimappower = NULL;
  } else
  {
    multimappower = bcktab_multimappower(dcov->bcktab);
  }
  esr = newEncodedsequencescanstate();
  diffptr = dcov->diffvalues;
  afterend = dcov->diffvalues + dcov->size;
  for (pos = 0; pos <= dcov->totallength; pos++)
  {
    gt_assert(modvalue == MODV(pos));
    gt_assert(diffptr == afterend || *diffptr >= (Diffvalue) modvalue);
    if (diffptr < afterend && (Diffvalue) modvalue == *diffptr)
    {
      if (withcheck)
      {
        idx = differencecover_packsamplepos(dcov,pos);
        gt_assert(idx < maxsamplesize && sampleidxused != NULL);
        if (ISIBITSET(sampleidxused,idx))
        {
          fprintf(stderr,"sample index %lu for pos %lu already used before\n",
                         idx,(unsigned long) pos);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        SETIBIT(sampleidxused,idx);
      }
      /* printf("pos mod %u in difference cover\n",dcov->vparam); */
      code = extractprefixcode(&unitsnotspecial,
                               dcov->encseq,
                               filltable,
                               readmode,
                               esr,
                               multimappower,
                               pos,
                               dcov->prefixlength);
      diffptr++;
    }
    if (modvalue < dcov->vmodmask)
    {
      modvalue++;
    } else
    {
      modvalue = 0;
      diffptr = dcov->diffvalues;
    }
  }
  printf("%lu positions are sampled (%.2f) wasted=%lu, pl=%u\n",
                                  dcov->samplesize,
                                  100.0 *
                                  (double) dcov->samplesize/dcov->totallength,
                                  maxsamplesize - dcov->samplesize,
                                  dcov->prefixlength);

  gt_free(filltable);
  if (esr != NULL)
  {
    freeEncodedsequencescanstate(&esr);
  }
  gt_free(sampleidxused);
}

void differencecovers_check(Seqpos maxcheck,const Encodedsequence *encseq,
                            Readmode readmode)
{
  Differencecover *dcov;
  size_t logmod, next = 0;
  unsigned int j, vparam;
  Seqpos pos1, pos2;

  if (maxcheck > getencseqtotallength(encseq))
  {
    maxcheck = getencseqtotallength(encseq);
  }
  for (logmod = 0;
       logmod < sizeof (differencecoversizes)/sizeof (differencecoversizes[0]);
       logmod++)
  {
    vparam = 1U << logmod;
    dcov = differencecover_new(vparam,encseq);
    if (dcov == NULL)
    {
      fprintf(stderr,"no difference cover for v=%u\n",vparam);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    for (j = 0; j<dcov->size; j++)
    {
      gt_assert(dcov->diffvalues[j] == differencecovertab[next]);
      next++;
    }
    for (pos1=0; pos1<maxcheck; pos1++)
    {
      for (pos2=0; pos2<maxcheck; pos2++)
      {
        gt_assert(differencecover_offset(dcov,pos1,pos2) ==
                  differencecover_offset(dcov,pos1,pos2));
      }
    }
    printf("v=%u (size=%u): ",dcov->vparam,dcov->size);
    differencecover_sample(dcov,readmode,true);
    differencecover_delete(dcov);
  }
  printf("# %u difference covers checked\n",(unsigned int) logmod);
  gt_assert(next == sizeof (differencecovertab)/sizeof (differencecovertab[0]));
}
