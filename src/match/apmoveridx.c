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
#include "core/symboldef.h"
#include "core/unused_api.h"
#include "core/ma.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "defined-types.h"
#include "apmoveridx.h"
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
#ifdef SKDEBUG
  unsigned long lastscorevalue;   /* the score for the given depth */
#endif
} Myerscolumn;

typedef struct
{
  bool skpp;
  unsigned long maxintervalwidth,
                patternlength,
                maxdistance,
                *eqsvector;
} Matchtaskinfo;

#ifdef SKDEBUG

static void showmaxleqvalue(FILE *fp,unsigned long maxleqk,
                            const Matchtaskinfo *mti)
{
  if (maxleqk == UNDEFMAXLEQK)
  {
    fprintf(fp,"undefined");
  } else
  {
    fprintf(fp,"%lu",maxleqk);
  }
}

static void apm_showMyerscolumn(const DECLAREPTRDFSSTATE(aliascol),
                                unsigned long score,
                                const void *dfsconstinfo)
{
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  const Myerscolumn *col = (const Myerscolumn *) aliascol;

  if (col->maxleqk == UNDEFMAXLEQK)
  {
    printf("[]");
  } else
  {
    unsigned long idx, backmask;

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

static void verifycolumnvalues(const Matchtaskinfo *mti,
                               const Myerscolumn *col,
                               unsigned long startscore)
{
  unsigned long idx, score, minscore, mask, bfmaxleqk;

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
    exit(EXIT_FAILURE);
  }
  if (bfmaxleqk != UNDEFMAXLEQK && minscore != col->lastscorevalue)
  {
    fprintf(stderr,"correct score = %lu != %lu = col->score\n",
                 minscore,
                 col->lastscorevalue);
    exit(EXIT_FAILURE);
  }
}

#endif

static void apm_initdfsconstinfo(void *dfsconstinfo,
                                 unsigned int alphasize,
                                 const Uchar *pattern,
                                 unsigned long patternlength,
                                 unsigned long maxdistance,
                                 unsigned long maxintervalwidth,
                                 bool skpp)
{
  Matchtaskinfo *mti = (Matchtaskinfo *) dfsconstinfo;

  initeqsvector(mti->eqsvector,(unsigned long) alphasize,pattern,patternlength);
  mti->patternlength = patternlength;
  mti->maxdistance = maxdistance;
  mti->maxintervalwidth = maxintervalwidth;
  mti->skpp = skpp;
}

static void *apm_allocatedfsconstinfo(unsigned int alphasize)
{
  Matchtaskinfo *mti = gt_malloc(sizeof(Matchtaskinfo));
  mti->eqsvector = gt_malloc(sizeof(*mti->eqsvector) * alphasize);
  return mti;
}

static void apm_freedfsconstinfo(void **dfsconstinfo)
{
  Matchtaskinfo *mti = (Matchtaskinfo *) *dfsconstinfo;

  gt_free(mti->eqsvector);
  gt_free(mti);
  *dfsconstinfo = NULL;
}

static void apm_initMyerscolumn(DECLAREPTRDFSSTATE(aliascolumn),
                                void *dfsconstinfo)
{
  Myerscolumn *column = (Myerscolumn *) aliascolumn;
  const Matchtaskinfo *mti = (Matchtaskinfo *) dfsconstinfo;

  column->Mv = 0UL;
  if (mti->skpp)
  {
    column->Pv = 0UL; /* first column consists of 0 => skip pattern prefix */
    column->maxleqk = mti->patternlength;
#ifdef SKDEBUG
    column->lastscorevalue = 0;
#endif
  } else
  {
    column->Pv = ~0UL; /* first column: 0 1 2 ... m */
    column->maxleqk = mti->maxdistance;
#ifdef SKDEBUG
    column->lastscorevalue = mti->maxdistance;
#endif
  }
}

static unsigned long apm_nextstepfullmatches(DECLAREPTRDFSSTATE(aliascolumn),
                                             GT_UNUSED Seqpos leftbound,
                                             GT_UNUSED Seqpos rightbound,
                                             Seqpos width,
                                             GT_UNUSED
                                             unsigned long currentdepth,
                                             void *dfsconstinfo)
{
  const Matchtaskinfo *mti = (Matchtaskinfo *) dfsconstinfo;
  Myerscolumn *col = (Myerscolumn *) aliascolumn;

  if (col->maxleqk == UNDEFMAXLEQK)
  {
    return 0; /* stop depth first traversal */
  }
  if (mti->maxintervalwidth == 0 || width == (Seqpos) 1)
  {
    if (col->maxleqk == SUCCESSMAXLEQK)
    {
      return mti->patternlength+1; /* success with match of length plen */
    }
  } else
  {
    if (width <= (Seqpos) mti->maxintervalwidth)
    {
      assert(col->maxleqk > 0);
      return col->maxleqk+1; /* success with match of length maxleqk */
    }
  }
  return 1UL; /* continue with depth first traversal */
}

/* implement types Limdfsstate and all function appearing until here */

static void apm_nextMyercolumn(const void *dfsconstinfo,
                               DECLAREPTRDFSSTATE(aliasoutcol),
                               GT_UNUSED unsigned long currentdepth,
                               Uchar currentchar,
                               const DECLAREPTRDFSSTATE(aliasincol))
{
  unsigned long Eq = 0, Xv, Xh, Ph, Mh, /* as in Myers Paper */
                backmask,               /* only one bit is on */
                idx,                    /* a counter */
                score;                  /* current score */
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  Myerscolumn *outcol = (Myerscolumn *) aliasoutcol;
  const Myerscolumn *incol = (const Myerscolumn *) aliasincol;

  assert(incol->maxleqk != UNDEFMAXLEQK);
  assert(mti->maxintervalwidth > 0 || incol->maxleqk != SUCCESSMAXLEQK);
  assert(currentchar != (Uchar) SEPARATOR);
  if (currentchar != (Uchar) WILDCARD)
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
#ifdef SKDEBUG
    outcol->lastscorevalue = incol->lastscorevalue;
#endif
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
#ifdef SKDEBUG
              outcol->lastscorevalue = score;
#endif
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
#ifdef SKDEBUG
      outcol->lastscorevalue = incol->lastscorevalue;
#endif
    }
  }
#ifdef SKDEBUG
  verifycolumnvalues(mti,outcol,currentdepth);
#endif
}

static void apm_inplacenextMyercolumn(const void *dfsconstinfo,
                                      DECLAREPTRDFSSTATE(aliascol),
                                      GT_UNUSED unsigned long currentdepth,
                                      Uchar currentchar)
{
  unsigned long Eq = 0, Xv, Xh, Ph, Mh, /* as in Myers Paper */
                backmask,           /* only one bit is on */
                idx,                /* a counter */
                score;              /* current score */
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  Myerscolumn *col = (Myerscolumn *) aliascol;

  assert(col->maxleqk != UNDEFMAXLEQK);
  assert(mti->maxintervalwidth > 0 || col->maxleqk != SUCCESSMAXLEQK);
  if (currentchar != (Uchar) WILDCARD)
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
#ifdef SKDEBUG
              col->lastscorevalue = score;
#endif
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

const AbstractDfstransformer *apm_AbstractDfstransformer(void)
{
  static const AbstractDfstransformer apm_adfst =
  {
    sizeof (Myerscolumn),
    apm_allocatedfsconstinfo,
    apm_initdfsconstinfo,
    NULL, /* no extractdfsconstinfo */
    apm_freedfsconstinfo,
    apm_initMyerscolumn,
    apm_nextstepfullmatches,
    apm_nextMyercolumn,
    apm_inplacenextMyercolumn,
#ifdef SKDEBUG
    apm_showMyerscolumn,
#endif
  };
  return &apm_adfst;
}
