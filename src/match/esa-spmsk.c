/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include "core/unused_api.h"
#include "core/readmode_api.h"
#include "core/error_api.h"
#include "core/arraydef.h"
#include "esa-bottomup.h"
#include "esa-spmsk.h"

typedef struct
{
  unsigned long firstinW;
} GtBUinfo_spmsk;

struct GtBUstate_spmsk /* global information */
{
  const GtEncseq *encseq;
  GtReadmode readmode;
  unsigned long totallength,
                minmatchlength,
                spmcounter;
  bool countspms,
       outputspms;
  GtArrayGtUlong Wset, Lset;
  /* Declare the stack as void as the real type
     GtArrayGtBUItvinfo_spmsk is declared later in esa-bottomup-spmsk.inc */
  void *stack;
};

static void initBUinfo_spmsk(GtBUinfo_spmsk *buinfo,
                             GT_UNUSED GtBUstate_spmsk *state)
{
  buinfo->firstinW = ULONG_MAX;
}

static void freeBUinfo_spmsk(GT_UNUSED GtBUinfo_spmsk *info,
                             GT_UNUSED GtBUstate_spmsk *state)
{
  return;
}

static int processleafedge_spmsk(bool firstedge,
                                 unsigned long fd,
                                 GT_UNUSED unsigned long flb,
                                 GtBUinfo_spmsk *finfo,
                                 unsigned long position,
                                 unsigned long seqnum,
                                 unsigned long relpos,
                                 GtBUstate_spmsk *state,
                                 GT_UNUSED GtError *err)

{
  if (fd >= state->minmatchlength)
  {
    unsigned long idx;

    if (firstedge)
    {
      gt_assert(finfo != NULL);
      ((GtBUinfo_spmsk *) finfo)->firstinW = state->Wset.nextfreeGtUlong;
    }
    if (position == 0 || gt_encseq_position_is_separator(state->encseq,
                                                         position - 1,
                                                         state->readmode))
    {
      gt_assert(relpos == 0);
      idx = gt_encseq_seqnum(state->encseq,position);
      GT_STOREINARRAY(&state->Wset,GtUlong,128,idx);
    } else
    {
      gt_assert(relpos > 0);
    }
    if (position + fd == state->totallength ||
        gt_encseq_position_is_separator(state->encseq,
                                        position + fd,state->readmode))
    {
      gt_assert(relpos + fd == gt_encseq_seqlength(state->encseq,seqnum));
      idx = gt_encseq_seqnum(state->encseq,position);
      gt_assert(idx == seqnum);
      GT_STOREINARRAY(&state->Lset,GtUlong,128,idx);
    } else
    {
      gt_assert(relpos + fd < gt_encseq_seqlength(state->encseq,seqnum));
    }
  }
  return 0;
}

static int processlcpinterval_spmsk(unsigned long lcp,
                                    GT_UNUSED unsigned long lb,
                                    GT_UNUSED unsigned long rb,
                                    GtBUinfo_spmsk *info,
                                    GtBUstate_spmsk *state,
                                    GT_UNUSED GtError *err)
{
  if (lcp >= state->minmatchlength)
  {
    unsigned long lidx, widx, firstpos;

    gt_assert(info != NULL);
    firstpos = ((GtBUinfo_spmsk *) info)->firstinW;
    for (lidx = 0; lidx < state->Lset.nextfreeGtUlong; lidx++)
    {
      if (state->outputspms)
      {
        unsigned long lpos = state->Lset.spaceGtUlong[lidx];

        for (widx = firstpos; widx < state->Wset.nextfreeGtUlong; widx++)
        {
          printf("%lu %lu %lu\n",lpos,state->Wset.spaceGtUlong[widx],lcp);
        }
      } else
      {
        gt_assert(state->countspms);
        if (firstpos < state->Wset.nextfreeGtUlong)
        {
          state->spmcounter += state->Wset.nextfreeGtUlong - firstpos;
        }
      }
    }
    state->Lset.nextfreeGtUlong = 0;
  } else
  {
    state->Wset.nextfreeGtUlong = 0;
  }
  return 0;
}

#define GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
#define GT_ESA_BOTTOM_UP_RAM
#include "esa-bottomup-spmsk.inc"

GtBUstate_spmsk *gt_spmsk_inl_new(const GtEncseq *encseq,
                            GtReadmode readmode,
                            unsigned long minmatchlength,
                            bool countspms,
                            bool outputspms,
                            GT_UNUSED const char *indexname)
{
  GtBUstate_spmsk *state = gt_malloc(sizeof (*state));

  state->encseq = encseq;
  state->readmode = readmode;
  state->totallength = gt_encseq_total_length(encseq);
  state->minmatchlength = minmatchlength;
  state->countspms = countspms;
  state->outputspms = outputspms;
  state->spmcounter = 0;
  state->stack = (void *) gt_GtArrayGtBUItvinfo_new();
  GT_INITARRAY(&state->Wset,GtUlong);
  GT_INITARRAY(&state->Lset,GtUlong);
  return state;
}

void gt_spmsk_inl_delete(GtBUstate_spmsk *state)
{
  if (state != NULL)
  {
    if (state->countspms)
    {
      printf("number of suffix-prefix matches=%lu\n",state->spmcounter);
    }
    GT_FREEARRAY(&state->Wset,GtUlong);
    GT_FREEARRAY(&state->Lset,GtUlong);
    gt_GtArrayGtBUItvinfo_delete_spmsk(
                 (GtArrayGtBUItvinfo_spmsk *) state->stack,state);
    gt_free(state);
  }
}

int gt_spmsk_inl_process(void *data,
                         const unsigned long *suftab_bucket,
                         const GtSeqnumrelpostab *snrp,
                         const uint16_t *lcptab_bucket,
                         unsigned long nonspecials,
                         GtError *err)
{
  GtBUstate_spmsk *state = (GtBUstate_spmsk *) data;

  if (gt_esa_bottomup_RAM_spmsk(suftab_bucket,
                                snrp,
                                lcptab_bucket,
                                nonspecials,
                                (GtArrayGtBUItvinfo_spmsk *) state->stack,
                                state,
                                err) != 0)
  {
    return -1;
  }
  return 0;
}
