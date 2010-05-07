/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include <stdarg.h>
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/ma_api.h"

#include "core/encseq.h"
#include "core/defined-types.h"
#include "apmeoveridx.h"
#include "absdfstrans-imp.h"
#include "initeqsvec.h"

#define UNDEFMAXLEQK      (mti->patternlength+1)
#define SUCCESSMAXLEQK    mti->patternlength

typedef struct
{
  unsigned long Pv,       /* the plus-vector for Myers Algorithm */
                Mv,       /* the minus-vector for Myers Algorithm */
                maxleqk;  /* \(\max\{i\in[0,m]\mid D(i)\leq k\}\) where
                          \(m\) is the length of the pattern, \(k\) is the
                          distance threshold, and \(D\) is
                          the current distance column */
  unsigned long lastdistvalue;   /* the score for the given depth */
} ApmeoveridxLimdfsstate;

typedef struct
{
  bool skpp;
  unsigned long maxintervalwidth,
                patternlength,
                maxdistance,
                *eqsvector;
} ApmeoveridxLimdfsconstinfo;

#ifdef SKDEBUG

static void showmaxleqvalue(FILE *fp,unsigned long maxleqk,
                            const Limdfsconstinfo *mti)
{
  if (maxleqk == UNDEFMAXLEQK)
  {
    fprintf(fp,"undefined");
  } else
  {
    fprintf(fp,"%lu",maxleqk);
  }
}

static void apme_showLimdfsstate(const DECLAREPTRDFSSTATE(aliascol),
                                 unsigned long currentdepth,
                                 const Limdfsconstinfo *mt)
{
  const ApmeoveridxLimdfsstate *col = (const ApmeoveridxLimdfsstate *) aliascol;
  const ApmeoveridxLimdfsconstinfo *mti =
                                        (const ApmeoveridxLimdfsconstinfo *) mt;

  if (col->maxleqk == UNDEFMAXLEQK)
  {
    printf("[]");
  } else
  {
    unsigned long idx, backmask, score = currentdepth;

    printf("[%lu",score);
    for (idx=1UL, backmask = 1UL; idx<=col->maxleqk; idx++, backmask <<= 1)
    {
      if (col->Pv & backmask)
      {
        score++;
      } else
      {
        if (col->Mv & backmask)
        {
          score--;
        }
      }
      printf(",%lu",score);
    }
    printf("] with maxleqk=%lu",col->maxleqk);
  }
}

static void verifycolumnvalues(const Limdfsconstinfo *mt,
                               const Limdfsstate *acol,
                               unsigned long startscore)
{
  unsigned long idx, score, minscore, mask, bfmaxleqk;
  const ApmeoveridxLimdfsstate *col = (const ApmeoveridxLimdfsstate *) acol;
  const ApmeoveridxLimdfsconstinfo *mti =
                                        (const ApmeoveridxLimdfsconstinfo *) mt;

  if (startscore <= mti->maxdistance)
  {
    bfmaxleqk = 0;
    minscore = startscore;
  } else
  {
    bfmaxleqk = UNDEFMAXLEQK;
    minscore = 0;
  }
  score = startscore;
  for (idx=1UL, mask = 1UL; idx <= mti->patternlength; idx++, mask <<= 1)
  {
    if (col->Pv & mask)
    {
      score++;
    } else
    {
      if (col->Mv & mask)
      {
        score--;
      }
    }
    if (score <= mti->maxdistance)
    {
      bfmaxleqk = idx;
      minscore = score;
    }
  }
  if (bfmaxleqk != col->maxleqk)
  {
    fprintf(stderr,"correct maxleqk = ");
    showmaxleqvalue(stderr,bfmaxleqk,mti);
    fprintf(stderr," != ");
    showmaxleqvalue(stderr,col->maxleqk,mti);
    fprintf(stderr," = col->maxleqk\n");
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (bfmaxleqk != UNDEFMAXLEQK && minscore != col->lastdistvalue)
  {
    fprintf(stderr,"correct score = %lu != %lu = col->score\n",
                 minscore,
                 col->lastdistvalue);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

#endif

static void apme_initdfsconstinfo(Limdfsconstinfo *mt,
                                  unsigned int alphasize,
                                  ...)
                                 /* Variable argument list is as follows:
                                    unsigned int alphasize
                                    const GtUchar *pattern,
                                    unsigned long patternlength,
                                    unsigned long maxdistance,
                                    unsigned long maxintervalwidth,
                                    bool skpp
                                 */
{
  va_list ap;
  const GtUchar *pattern;
  ApmeoveridxLimdfsconstinfo *mti = (ApmeoveridxLimdfsconstinfo *) mt;

  va_start(ap,alphasize);
  pattern = va_arg(ap, const GtUchar *);
  mti->patternlength = va_arg(ap, unsigned long);
  mti->maxdistance = va_arg(ap, unsigned long);
  mti->maxintervalwidth = va_arg(ap, unsigned long);
  mti->skpp = (bool) va_arg(ap, int);
  va_end(ap);
  gt_assert(mti->maxdistance < mti->patternlength);
  gt_initeqsvector(mti->eqsvector,(unsigned long) alphasize,
                pattern,mti->patternlength);
}

static Limdfsconstinfo *apme_allocatedfsconstinfo(unsigned int alphasize)
{
  ApmeoveridxLimdfsconstinfo *mti =
                                 gt_malloc(sizeof (ApmeoveridxLimdfsconstinfo));
  mti->eqsvector = gt_malloc(sizeof (*mti->eqsvector) * alphasize);
  return (Limdfsconstinfo*) mti;
}

static void apme_freedfsconstinfo(Limdfsconstinfo **mtptr)
{
  ApmeoveridxLimdfsconstinfo **mtiptr = (ApmeoveridxLimdfsconstinfo **) mtptr;
  gt_free((*mtiptr)->eqsvector);
  gt_free(*mtiptr);
  *mtiptr = NULL;
}

static void apme_initLimdfsstate(DECLAREPTRDFSSTATE(aliascolumn),
                                 Limdfsconstinfo *mt)
{
  ApmeoveridxLimdfsstate *column = (ApmeoveridxLimdfsstate *) aliascolumn;
  ApmeoveridxLimdfsconstinfo *mti = (ApmeoveridxLimdfsconstinfo *) mt;

  column->Mv = 0UL;
  if (mti->skpp)
  {
    column->Pv = 0UL; /* first column consists of 0 => skip pattern prefix */
    column->maxleqk = mti->patternlength;
    column->lastdistvalue = 0;
  } else
  {
    column->Pv = ~0UL; /* first column: 0 1 2 ... m */
    column->maxleqk = mti->maxdistance;
    column->lastdistvalue = mti->maxdistance;
  }
}

static void apme_fullmatchLimdfsstate(Limdfsresult *limdfsresult,
                                      DECLAREPTRDFSSTATE(aliascolumn),
                                      GT_UNUSED unsigned long leftbound,
                                      GT_UNUSED unsigned long rightbound,
                                      unsigned long width,
                                      GT_UNUSED unsigned long currentdepth,
                                      Limdfsconstinfo *mt)
{
  ApmeoveridxLimdfsstate *col = (ApmeoveridxLimdfsstate *) aliascolumn;
  ApmeoveridxLimdfsconstinfo *mti = (ApmeoveridxLimdfsconstinfo *) mt;

  if (col->maxleqk == UNDEFMAXLEQK)
  {
    limdfsresult->status = Limdfsstop; /* stop depth first traversal */
    return;
  }
  if (mti->maxintervalwidth == 0 || width == (unsigned long) 1)
  {
    if (col->maxleqk == SUCCESSMAXLEQK)
    {
      /* success with match of length plen */
      limdfsresult->status = Limdfssuccess;
      limdfsresult->pprefixlen = mti->patternlength;
      limdfsresult->distance = col->lastdistvalue;
      return;
    }
  } else
  {
    if (width <= (unsigned long) mti->maxintervalwidth)
    {
      /* success with match of length maxleqk */
      gt_assert(col->maxleqk > 0);
      limdfsresult->status = Limdfssuccess;
      limdfsresult->pprefixlen = col->maxleqk;
      limdfsresult->distance = col->lastdistvalue;
      return;
    }
  }
  /* continue with depth first traversal */
  limdfsresult->status = Limdfscontinue;
}

static void apme_nextLimdfsstate(const Limdfsconstinfo *mt,
                                 DECLAREPTRDFSSTATE(aliasoutcol),
                                 GT_UNUSED unsigned long currentdepth,
                                 GtUchar currentchar,
                                 const DECLAREPTRDFSSTATE(aliasincol))
{
  unsigned long Eq = 0, Xv, Xh, Ph, Mh, /* as in Myers Paper */
                backmask,               /* only one bit is on */
                idx,                    /* a counter */
                score;                  /* current score */
  ApmeoveridxLimdfsstate *outcol = (ApmeoveridxLimdfsstate *) aliasoutcol;
  const ApmeoveridxLimdfsstate *incol =
                                    (const ApmeoveridxLimdfsstate *) aliasincol;
  const ApmeoveridxLimdfsconstinfo *mti =
                                        (const ApmeoveridxLimdfsconstinfo *) mt;

  gt_assert(incol->maxleqk != UNDEFMAXLEQK);
  gt_assert(mti->maxintervalwidth > 0 || incol->maxleqk != SUCCESSMAXLEQK);
  gt_assert(currentchar != (GtUchar) SEPARATOR);
  if (currentchar != (GtUchar) WILDCARD)
  {
    Eq = mti->eqsvector[(unsigned long) currentchar];
  }
  Xv = Eq | incol->Mv;
  Xh = (((Eq & incol->Pv) + incol->Pv) ^ incol->Pv) | Eq;

  Ph = incol->Mv | ~ (Xh | incol->Pv);
  Mh = incol->Pv & Xh;

  Ph = (Ph << 1) | 1UL;
  outcol->Pv = (Mh << 1) | ~ (Xv | Ph);
  outcol->Mv = Ph & Xv;
  backmask = 1UL << incol->maxleqk;
  if (Eq & backmask || Mh & backmask)
  {
    outcol->maxleqk = incol->maxleqk + 1UL;
    outcol->lastdistvalue = incol->lastdistvalue;
  } else
  {
    if (Ph & backmask)
    {
      score = mti->maxdistance+1;
      outcol->maxleqk = UNDEFMAXLEQK;
      if (incol->maxleqk > 0)
      {
        for (idx = incol->maxleqk - 1, backmask >>= 1;
             /* Nothing */;
             backmask >>= 1)
        {
          if (outcol->Pv & backmask)
          {
            score--;
            if (score <= mti->maxdistance)
            {
              outcol->maxleqk = idx;
              outcol->lastdistvalue = score;
              break;
            }
          } else
          {
            if (outcol->Mv & backmask)
            {
              score++;
            }
          }
          if (idx > 0)
          {
            idx--;
          } else
          {
            break;
          }
        }
      }
    } else
    {
      outcol->maxleqk = incol->maxleqk;
      outcol->lastdistvalue = incol->lastdistvalue;
    }
  }
#ifdef SKDEBUG
  verifycolumnvalues(mti,outcol,currentdepth);
#endif
}

static void apme_inplacenextLimdfsstate(const Limdfsconstinfo *mt,
                                        DECLAREPTRDFSSTATE(aliascol),
                                        GT_UNUSED unsigned long currentdepth,
                                        GtUchar currentchar)
{
  unsigned long Eq = 0, Xv, Xh, Ph, Mh, /* as in Myers Paper */
                backmask,           /* only one bit is on */
                idx,                /* a counter */
                score;              /* current score */
  ApmeoveridxLimdfsstate *col = (ApmeoveridxLimdfsstate *) aliascol;
  const ApmeoveridxLimdfsconstinfo *mti =
                                        (const ApmeoveridxLimdfsconstinfo *) mt;

  gt_assert(col->maxleqk != UNDEFMAXLEQK);
  gt_assert(mti->maxintervalwidth > 0 || col->maxleqk != SUCCESSMAXLEQK);
  if (currentchar != (GtUchar) WILDCARD)
  {
    Eq = mti->eqsvector[(unsigned long) currentchar];
  }
  Xv = Eq | col->Mv;
  Xh = (((Eq & col->Pv) + col->Pv) ^ col->Pv) | Eq;

  Ph = col->Mv | ~ (Xh | col->Pv);
  Mh = col->Pv & Xh;

  Ph = (Ph << 1) | 1UL;
  col->Pv = (Mh << 1) | ~ (Xv | Ph);
  col->Mv = Ph & Xv;
  backmask = 1UL << col->maxleqk;
  if (Eq & backmask || Mh & backmask)
  {
    col->maxleqk++;
  } else
  {
    if (Ph & backmask)
    {
      unsigned long tmpmaxleqk = UNDEFMAXLEQK;
      score = mti->maxdistance+1;
      if (col->maxleqk > 0)
      {
        for (idx = col->maxleqk - 1, backmask >>= 1;
             /* Nothing */;
             backmask >>= 1)
        {
          if (col->Pv & backmask)
          {
            score--;
            if (score <= mti->maxdistance)
            {
              tmpmaxleqk = idx;
              col->lastdistvalue = score;
              break;
            }
          } else
          {
            if (col->Mv & backmask)
            {
              score++;
            }
          }
          if (idx > 0)
          {
            idx--;
          } else
          {
            break;
          }
        }
      }
      col->maxleqk = tmpmaxleqk;
    }
  }
}

const AbstractDfstransformer *gt_apme_AbstractDfstransformer(void)
{
  static const AbstractDfstransformer apme_adfst =
  {
    sizeof (ApmeoveridxLimdfsstate),
    apme_allocatedfsconstinfo,
    apme_initdfsconstinfo,
    NULL, /* no extractdfsconstinfo */
    apme_freedfsconstinfo,
    apme_initLimdfsstate,
    NULL,
    NULL,
    NULL,
    apme_fullmatchLimdfsstate,
    apme_nextLimdfsstate,
    apme_inplacenextLimdfsstate,
#ifdef SKDEBUG
    apme_showLimdfsstate,
#endif
  };
  return &apme_adfst;
}
