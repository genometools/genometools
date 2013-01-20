/*
  Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh-uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#include "core/progressbar.h"
#include "match/esa-seqread.h"
#include "match/esa-bottomup.h"
#include "match/rdj-revcompl-def.h"
#include "match/rdj-contfind-bottomup.h"

typedef struct
{
  /* length of the shortest sequence in the array */
  unsigned long shortest;

  unsigned long spacing; /* for eqlen read length mode */

  const GtEncseq *encseq;
  GtBitsequence *sspbittab;

  GtBitsequence *contained;

  unsigned long cmin;
  unsigned long csize;

  /* progressbar */
  bool show_progressbar;
  unsigned long long progress;

  /* revcompl mode */
  unsigned long firstrevcompl;
  unsigned long nofsequences;

  unsigned long counter;
} ContfindBUstate;

typedef ContfindBUstate GtBUstate_rdjce;
typedef ContfindBUstate GtBUstate_rdjcv;
typedef struct {} GtBUinfo_rdjce;
typedef struct {} GtBUinfo_rdjcv;

static void initBUinfo_rdjcv(GT_UNUSED GtBUinfo_rdjcv *info,
    GT_UNUSED GtBUstate_rdjcv *state)
{
  /* nothing to do */
}

static void initBUinfo_rdjce(GT_UNUSED GtBUinfo_rdjce *info,
    GT_UNUSED GtBUstate_rdjce *state)
{
  /* nothing to do */
}

static inline void freeBUinfo_rdjcv(GT_UNUSED GtBUinfo_rdjcv *info,
    GT_UNUSED GtBUstate_rdjcv *state)
{
  /* nothing to do */
}

static inline void freeBUinfo_rdjce(GT_UNUSED GtBUinfo_rdjce *info,
    GT_UNUSED GtBUstate_rdjce *state)
{
  /* nothing to do */
}

static inline void processcontained(unsigned long seqnum,
    ContfindBUstate *state)
{
  if (state->firstrevcompl > 0)
    seqnum = GT_READJOINER_READNUM(seqnum, state->firstrevcompl,
        state->nofsequences);
  if (!GT_ISIBITSET(state->contained, seqnum))
  {
    GT_SETIBIT(state->contained, seqnum);
    state->counter++;
  }
  if (state->csize == 0 || seqnum < state->cmin)
    state->cmin = seqnum;
  state->csize++;
}

static inline int processleafedge_rdjcv(GT_UNUSED bool firstsucc,
    unsigned long fatherdepth,
    GT_UNUSED GtBUinfo_rdjcv *father, unsigned long leafnumber,
    GtBUstate_rdjcv *state, GT_UNUSED GtError *err)
{
  unsigned long seqnum;

  if (fatherdepth >= state->shortest)
  {
    if ((leafnumber == 0 ||
        GT_ISIBITSET(state->sspbittab, leafnumber-1)) &&
        GT_ISIBITSET(state->sspbittab, leafnumber + fatherdepth))
    {
      seqnum = gt_encseq_seqnum(state->encseq, leafnumber);
      processcontained(seqnum, state);
    }
  }
  if (state->show_progressbar) state->progress++;
  return 0;
}

static inline int processleafedge_rdjce(GT_UNUSED bool firstsucc,
    unsigned long fatherdepth,
    GT_UNUSED GtBUinfo_rdjce *father, unsigned long leafnumber,
    GtBUstate_rdjce *state, GT_UNUSED GtError *err)
{
  unsigned long seqnum;

  if (fatherdepth == state->shortest && (leafnumber % state->spacing) == 0)
  {
    seqnum = leafnumber / state->spacing;
    processcontained(seqnum, state);
  }
  if (state->show_progressbar) state->progress++;
  return 0;
}

static inline int processbranchingedge_rdjce(GT_UNUSED bool firstsucc,
    GT_UNUSED unsigned long fatherdepth,
    GT_UNUSED GtBUinfo_rdjce *father, unsigned long sondepth,
    GT_UNUSED unsigned long sonwidth,
    GT_UNUSED GtBUinfo_rdjce *son, GtBUstate_rdjce *state,
    GT_UNUSED GtError *err)
{
  if (sondepth == state->shortest /* == read length */)
  {
    /* keep one copy */
    if (state->csize > 0)
    {
      GT_UNSETIBIT(state->contained, state->cmin);
      state->counter--;
    }
  }
  state->csize = 0;
  return 0;
}

static inline int processbranchingedge_rdjcv(GT_UNUSED bool firstsucc,
    GT_UNUSED unsigned long fatherdepth,
    GT_UNUSED GtBUinfo_rdjcv *father, unsigned long sondepth,
    unsigned long sonwidth, GT_UNUSED GtBUinfo_rdjcv *son,
    GtBUstate_rdjcv *state, GT_UNUSED GtError *err)
{
  if (sondepth >= state->shortest)
  {
    /* keep one copy if several copies of a sequence are present
       and the sequence is not contained in any other */
    if (state->csize == sonwidth)
    {
      gt_assert(state->csize > 0);
      GT_UNSETIBIT(state->contained, state->cmin);
      state->counter--;
    }
  }
  state->csize = 0;
  return 0;
}

#include "match/esa-bottomup-rdjcv.inc"

#include "match/esa-bottomup-rdjce.inc"

/* prepare sspbittab and determine length of shortest sequence */
static void prepare_sspbittab_and_shortest(unsigned long totallength,
    ContfindBUstate *state)
{
  unsigned long length, lastseqstart, i, ssp;

  GT_INITBITTAB(state->sspbittab, totallength + 1);
  lastseqstart = 0;
  state->shortest = totallength;
  for (i = 1UL; i <= state->nofsequences - 1; i++)
  {
    ssp = gt_encseq_seqstartpos(state->encseq, i) - 1;
    GT_SETIBIT(state->sspbittab, ssp);
    length = ssp - lastseqstart;
    lastseqstart = ssp + 1;
    if (length < state->shortest)
      state->shortest = length;
  }
  GT_SETIBIT(state->sspbittab, totallength);
  length = totallength - lastseqstart;
  if (length < state->shortest)
    state->shortest = length;
}

unsigned long gt_contfind_bottomup(Sequentialsuffixarrayreader *ssar,
                     bool show_progressbar, GtBitsequence *contained,
                     unsigned long firstrevcompl,
                     unsigned long read_length /* 0 = variable */)
{
  ContfindBUstate state;
  unsigned long totallength;
  GT_UNUSED int retval;

  gt_assert(ssar != NULL);
  gt_assert(contained != NULL);

  state.contained = contained;
  state.encseq = gt_encseqSequentialsuffixarrayreader(ssar);
  totallength = gt_encseq_total_length(state.encseq);
  state.nofsequences = gt_encseq_num_of_sequences(state.encseq);

  if (read_length == 0)
  {
    prepare_sspbittab_and_shortest(totallength, &state);
  }
  else
  {
    state.shortest = read_length;
    state.spacing = read_length + 1;
  }

  state.show_progressbar = show_progressbar;
  state.csize            = 0;
  state.cmin             = 0;
  state.firstrevcompl    = firstrevcompl;
  state.counter          = 0;

  if (show_progressbar)
  {
    state.progress = 0;
    gt_progressbar_start(&(state.progress),
        (unsigned long long)totallength);
  }

  retval = (read_length == 0)
      ? gt_esa_bottomup_rdjcv(ssar, &state, NULL)
      : gt_esa_bottomup_rdjce(ssar, &state, NULL);
  gt_assert(retval == 0);

  if (show_progressbar)
    gt_progressbar_stop();
  if (read_length == 0)
    gt_free(state.sspbittab);
  return state.counter;
}
