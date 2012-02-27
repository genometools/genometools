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
                spmcounter,
                spaceforbucketprocessing;
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
                                 GtBUinfo_spmsk *finfo,
                                 unsigned long seqnum,
                                 unsigned long relpos,
                                 GtBUstate_spmsk *state,
                                 GT_UNUSED GtError *err)

{
  if (fd >= state->minmatchlength)
  {
    if (firstedge)
    {
      gt_assert(finfo != NULL);
      ((GtBUinfo_spmsk *) finfo)->firstinW = state->Wset.nextfreeGtUlong;
    }
    if (relpos == 0)
    {
      GT_STOREINARRAY(&state->Wset,GtUlong,128,seqnum);
    }
    if (relpos + fd == gt_encseq_seqlength(state->encseq,seqnum))
    {
      GT_STOREINARRAY(&state->Lset,GtUlong,128,seqnum);
    }
  }
  return 0;
}

static int processlcpinterval_spmsk(unsigned long lcp,
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

unsigned long gt_spmsk_inl_delete(GtBUstate_spmsk *state)
{
  if (state != NULL)
  {
    unsigned long tmpcount;

    GT_FREEARRAY(&state->Wset,GtUlong);
    GT_FREEARRAY(&state->Lset,GtUlong);
    gt_GtArrayGtBUItvinfo_delete_spmsk(
                 (GtArrayGtBUItvinfo_spmsk *) state->stack,state);
    tmpcount = state->spmcounter;
    gt_free(state);
    return tmpcount;
  }
  return 0;
}

int gt_spmsk_inl_process(void *data,
                         const unsigned long *seqnum_relpos_bucket,
                         const GtSeqnumrelpos *snrp,
                         const uint16_t *lcptab_bucket,
                         unsigned long nonspecials,
                         unsigned long spaceforbucketprocessing,
                         GtError *err)
{
  GtBUstate_spmsk *state = (GtBUstate_spmsk *) data;

  gt_assert(snrp != NULL);
  gt_assert(lcptab_bucket != NULL);
  state->spaceforbucketprocessing = spaceforbucketprocessing;
  if (gt_esa_bottomup_RAM_spmsk(seqnum_relpos_bucket,
                                lcptab_bucket,
                                nonspecials,
                                (GtArrayGtBUItvinfo_spmsk *) state->stack,
                                state,
                                snrp,
                                err) != 0)
  {
    return -1;
  }
  return 0;
}

void gt_spmsk_inl_process_end(GT_UNUSED void *data)
{
  return;
}
