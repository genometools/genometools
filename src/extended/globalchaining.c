/*
  Copyright (c) 2007, 2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004       Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2004, 2007 Center for Bioinformatics, University of Hamburg

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
#include "core/array.h"
#include "core/assert_api.h"
#include "core/ma.h"
#include "core/log.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "extended/globalchaining.h"

#define UNDEFPREVIOUS           num_of_fragments
#define GETSTOREDSTARTPOINT(DIM,IDX)\
        fragments[IDX].startpos##DIM
#define GETSTOREDENDPOINT(DIM,IDX)\
        fragments[IDX].endpos##DIM
#define GETSTOREDLENGTH(DIM,IDX)\
        (GETSTOREDENDPOINT(DIM,IDX) - GETSTOREDSTARTPOINT(DIM,IDX))

typedef struct {
  long maxscore;
  unsigned long maxfragnum;
  bool defined;
} Maxfragvalue;

typedef struct {
  unsigned long previousinchain;  /* previous index in chain */
  long score;                     /* score of highest scoring chain ending here,
                                     computed */
} GtChaininfo;

/*
  The following structure is used to store additional information for every
  fragment when chains containing overlaps are computed.
*/

typedef struct {
  bool active;                     /* fragment is active, i.e., it is an end
                                      fragment */
  unsigned long startofchain,      /* fragment number starting this chain */
                dim1lengthofchain, /* length of chain in the first dimension
                                      (necessary to compute the coverage) */
                chainarray;        /* contains computed chains:
                                      [startofchain] -> end fragment number */
} Overlapinfo;

static void initoverlapinfo(Overlapinfo *overlapinfo,
                            GtFragment *fragments,
                            unsigned long num_of_fragments)
{
  unsigned long i;

  for (i = 0; i < num_of_fragments; i++) {
    overlapinfo[i].active            = true;
    overlapinfo[i].startofchain      = i;
    overlapinfo[i].dim1lengthofchain = GETSTOREDLENGTH(1, i);
    overlapinfo[i].chainarray        = UNDEFPREVIOUS;
  }
}

static bool colinearfragments(GtFragment *fragments,
                              unsigned long i, unsigned long j)
{
  if (GETSTOREDSTARTPOINT(1, i) < GETSTOREDSTARTPOINT(1, j) &&
      GETSTOREDENDPOINT(1, i)   < GETSTOREDENDPOINT(1, j)   &&
      GETSTOREDSTARTPOINT(2, i) < GETSTOREDSTARTPOINT(2, j) &&
      GETSTOREDENDPOINT(2, i)   < GETSTOREDENDPOINT(2, j)) {
    return true;
  }
  return false;
}

static long overlapcost(GtFragment *fragments,
                        unsigned long i, unsigned long j)
{
  unsigned long overlaplength = 0;

  /* add overlap in first dimension */
  if (GETSTOREDSTARTPOINT(1, j) <= GETSTOREDENDPOINT(1, i))
    overlaplength += GETSTOREDENDPOINT(1, i) - GETSTOREDSTARTPOINT(1, j) + 1;

  /* add overlap in second dimension */
  if (GETSTOREDSTARTPOINT(2, j) <= GETSTOREDENDPOINT(2, i))
    overlaplength += GETSTOREDENDPOINT(2, i) - GETSTOREDSTARTPOINT(2, j) + 1;

  gt_log_log("overlap total  (#%lu, #%lu)=%lu", i, j, overlaplength);

  return (long) overlaplength;
}

static void chainingboundarycases(GtChain *chain,
                                  GtFragment *fragments,
                                  unsigned long num_of_fragments)
{
  if (num_of_fragments == 0)
    gt_chain_reset(chain);
  else if (num_of_fragments == 1UL) {
    gt_chain_reset(chain);
    gt_chain_set_score(chain, fragments[0].weight);
    gt_chain_add_fragnum(chain, 0);
  }
}

static void retracepreviousinchain(GtChain *chain, GtChaininfo *chaininfo,
                                   unsigned long num_of_fragments,
                                   unsigned long retracestart)
{
  unsigned long fragnum, idx, lengthofchain;

  for (lengthofchain = 0, fragnum = retracestart;
       fragnum != UNDEFPREVIOUS; lengthofchain++) {
    gt_chain_add_fragnum(chain, UNDEFPREVIOUS); /* add dummy */
    fragnum = chaininfo[fragnum].previousinchain;
  }
  fragnum = retracestart;
  idx = lengthofchain;
  while (fragnum != UNDEFPREVIOUS) {
    gt_assert(idx != 0);
    idx--;
    gt_chain_set_fragnum(chain, idx, fragnum); /* set dummy */
    fragnum = chaininfo[fragnum].previousinchain;
  }
  gt_assert(idx == 0);
  gt_assert(gt_chain_size(chain) == lengthofchain);
}

static bool check_max_gap_width(GtFragment *fragments,
                                unsigned long max_gap_width,
                                unsigned long leftfrag,
                                unsigned long rightfrag)
{
  unsigned long gapwidth, startpoint, endpoint;

  startpoint = GETSTOREDSTARTPOINT(1,rightfrag);
  endpoint = GETSTOREDENDPOINT(1,leftfrag);
  if (startpoint <= endpoint)
    gapwidth = 0;
  else
    gapwidth = startpoint - endpoint - 1;
  if (gapwidth > max_gap_width)
    return false;

  startpoint = GETSTOREDSTARTPOINT(2,rightfrag);
  endpoint = GETSTOREDENDPOINT(2,leftfrag);
  if (startpoint <= endpoint)
    gapwidth = 0;
  else
    gapwidth = startpoint - endpoint - 1;
  if (gapwidth > max_gap_width)
    return false;

  return true;
}

static void bruteforcechainingscores(GtChaininfo *chaininfo,
                                     unsigned long max_gap_width,
                                     GtFragment *fragments,
                                     unsigned long num_of_fragments,
                                     Overlapinfo *overlapinfo)
{
  unsigned long previous, leftfrag, rightfrag, overlaplength;
  long weightright, score;
  Maxfragvalue localmaxfrag;
  bool combinable;

  /* get rid of compiler warnings */
  localmaxfrag.maxscore = 0;
  localmaxfrag.maxfragnum = 0;

  if (num_of_fragments > 1UL) {
    chaininfo[0].previousinchain = UNDEFPREVIOUS;
    chaininfo[0].score = fragments[0].weight;
    for (rightfrag = 1UL; rightfrag < num_of_fragments; rightfrag++) {
      weightright = fragments[rightfrag].weight;
      localmaxfrag.defined = false;
      for (leftfrag = 0; leftfrag < rightfrag; leftfrag++) {
        if (max_gap_width != 0 &&
            !check_max_gap_width(fragments, max_gap_width, leftfrag,
                                 rightfrag)) {
          combinable = false;
        }
        else {
          combinable = colinearfragments(fragments, leftfrag, rightfrag);
        }
        if (combinable) {
          score = chaininfo[leftfrag].score
                  - overlapcost(fragments, leftfrag, rightfrag);
          if (score > 0) {
            score += weightright;
            previous = leftfrag;
          }
          else {
            score = weightright;
            previous = UNDEFPREVIOUS;
          }
          if (!localmaxfrag.defined || localmaxfrag.maxscore < score) {
            localmaxfrag.maxscore = score;
            localmaxfrag.maxfragnum = previous;
            localmaxfrag.defined = true;
          }
        }
      }
      if (localmaxfrag.defined) {
        chaininfo[rightfrag].previousinchain = localmaxfrag.maxfragnum;
        chaininfo[rightfrag].score = localmaxfrag.maxscore;
        if (overlapinfo && localmaxfrag.maxfragnum != UNDEFPREVIOUS) {
          overlapinfo[localmaxfrag.maxfragnum].active = false;
          overlapinfo[rightfrag].startofchain = overlapinfo[localmaxfrag
                                                            .maxfragnum]
                                                .startofchain;
          if (GETSTOREDSTARTPOINT(1, rightfrag) <=
              GETSTOREDENDPOINT(1, localmaxfrag.maxfragnum)) {
            overlaplength = GETSTOREDENDPOINT(1, localmaxfrag.maxfragnum)
                            - GETSTOREDSTARTPOINT(1, rightfrag) + 1;
          }
          else {
            overlaplength = 0;
          }
          overlapinfo[rightfrag].dim1lengthofchain += overlapinfo[localmaxfrag
                                                                  .maxfragnum]
                                                      .dim1lengthofchain
                                                      - overlaplength;
        }
      }
      else {
        chaininfo[rightfrag].previousinchain = UNDEFPREVIOUS;
        chaininfo[rightfrag].score = weightright;
      }
    }
  }
}

static bool isrightmaximallocalchain(GtChaininfo *chaininfo,
                                     unsigned long num_of_fragments,
                                     unsigned long currentfrag)
{
  if (currentfrag == num_of_fragments - 1)
    return true;
  if (chaininfo[currentfrag+1].previousinchain != currentfrag)
    return true;
  if (chaininfo[currentfrag+1].score < chaininfo[currentfrag].score)
    return true;
  return false;
}

static bool retrievemaximalscore(long *maxscore, GtChaininfo *chaininfo,
                                 unsigned long num_of_fragments)
{
  unsigned long i;
  bool maxscoredefined = false;

  *maxscore = 0;
  for (i=0; i< num_of_fragments; i++) {
    if (isrightmaximallocalchain(chaininfo, num_of_fragments,i)) {
      if (!maxscoredefined || *maxscore < chaininfo[i].score) {
        *maxscore = chaininfo[i].score;
        maxscoredefined = true;
      }
    }
  }
  return maxscoredefined;
}

static void retrievechainthreshold(GtChaininfo *chaininfo,
                                   GtFragment *fragments,
                                   unsigned long num_of_fragments,
                                   unsigned long max_gap_width,
                                   GtChain *chain, long minscore,
                                   GtChainProc chainprocessor,
                                   void *cpinfo)
{
  unsigned long i;

  for (i = 0; i < num_of_fragments; i++) {
    if (isrightmaximallocalchain(chaininfo, num_of_fragments,i)) {
      if (chaininfo[i].score >= minscore) {
        gt_chain_reset(chain);
        gt_chain_set_score(chain, chaininfo[i].score);
        retracepreviousinchain(chain, chaininfo, num_of_fragments, i);
        chainprocessor(chain, fragments, num_of_fragments, max_gap_width,
                       cpinfo);
      }
    }
  }
}

static void findmaximalscores(GtChain *chain, GtChaininfo *chaininfo,
                              GtFragment *fragments,
                              unsigned long num_of_fragments,
                              unsigned long max_gap_width,
                              GtChainProc chainprocessor, void *cpinfo)
{
  long minscore;
  bool minscoredefined = false;

  minscoredefined = retrievemaximalscore(&minscore, chaininfo,
                                         num_of_fragments);

  if (minscoredefined) {
    retrievechainthreshold(chaininfo, fragments, num_of_fragments,
                           max_gap_width, chain, minscore, chainprocessor,
                           cpinfo);
  }
}

static void findmaximalscores_withoverlaps(GtChain *chain,
                                           GtChaininfo *chaininfo,
                                           GtFragment *fragments,
                                           unsigned long num_of_fragments,
                                           unsigned long max_gap_width,
                                           unsigned long seqlen1,
                                           double mincoverage,
                                           GtChainProc chainprocessor,
                                           void *cpinfo,
                                           Overlapinfo *overlapinfo)
{
  unsigned long i, startfrag;
  GtArray *startfragments;

  gt_assert(seqlen1 != GT_UNDEF_ULONG);
  gt_assert(mincoverage != GT_UNDEF_DOUBLE);
  startfragments = gt_array_new(sizeof (unsigned long));

  /* compute chain array */
  for (i = 0; i < num_of_fragments; i++) {
    if (overlapinfo[i].active) {
      /* current fragment is active */
      if (overlapinfo[overlapinfo[i].startofchain].chainarray ==
          UNDEFPREVIOUS) {
        if (((double) overlapinfo[i].dim1lengthofchain / (double) seqlen1) >=
            mincoverage) {
          /* no other fragment has the same start fragment yet and coverage is
             high enough -> store end fragment */
          overlapinfo[overlapinfo[i].startofchain].chainarray = i;

          /* since this is the first time, store start fragment number
             to avoid additional scan of all fragments below */
          gt_array_add(startfragments, overlapinfo[i].startofchain);
        }
      }
      else if (overlapinfo[i].dim1lengthofchain >
               overlapinfo[overlapinfo[overlapinfo[i].startofchain].chainarray]
                          .dim1lengthofchain) {
        /* coverage is higher then coverage of earlier start fragment
           -> update end fragment */
        overlapinfo[overlapinfo[i].startofchain].chainarray = i;
      }
    }
  }

  /* retrieve maximal chains */
  for (i = 0; i < gt_array_size(startfragments); i++) {
    startfrag = *(unsigned long*) gt_array_get(startfragments, i);
    gt_assert(overlapinfo[startfrag].chainarray != UNDEFPREVIOUS);
    gt_chain_reset(chain);
    gt_chain_set_score(chain,
                       chaininfo[overlapinfo[startfrag].chainarray].score);
    retracepreviousinchain(chain, chaininfo, num_of_fragments,
                           overlapinfo[startfrag].chainarray);
    chainprocessor(chain, fragments, num_of_fragments, max_gap_width, cpinfo);
  }

  gt_array_delete(startfragments);
}

static void log_fragments(GtFragment *fragments, unsigned long num_of_fragments)
{
  unsigned long i;
  gt_log_log("show chaining fragments");
  for (i = 0; i < num_of_fragments; i++) {
    GtFragment *frag = fragments + i;
    gt_log_log("#%lu: s1=%lu, s1=%lu, l1=%lu, s2=%lu, e2=%lu, l2=%lu, w=%lu", i,
            frag->startpos1, frag->endpos1, frag->endpos1 - frag->startpos1 + 1,
            frag->startpos2, frag->endpos2, frag->endpos2 - frag->startpos2 + 1,
            frag->weight);
  }
}

static void globalchaining_generic(bool maxscore_chains,
                                   unsigned long max_gap_width,
                                   GtFragment *fragments,
                                   unsigned long num_of_fragments,
                                   unsigned long seqlen1, double mincoverage,
                                   GtChainProc chainprocessor, void *cpinfo)
{
  Overlapinfo *overlapinfo = NULL;
  GtChaininfo *chaininfo;
  GtChain *chain;
  chain = gt_chain_new();
  chaininfo = gt_malloc(sizeof (GtChaininfo) * num_of_fragments);
  if (gt_log_enabled())
    log_fragments(fragments, num_of_fragments);
  if (num_of_fragments > 1UL) {
    /* compute chains */
    if (!maxscore_chains) {
      overlapinfo = gt_malloc(sizeof (Overlapinfo) * num_of_fragments);
      initoverlapinfo(overlapinfo, fragments, num_of_fragments);
    }
    bruteforcechainingscores(chaininfo, max_gap_width, fragments,
                             num_of_fragments, overlapinfo);
    if (maxscore_chains) {
      findmaximalscores(chain, chaininfo, fragments, num_of_fragments,
                        max_gap_width, chainprocessor, cpinfo);
    }
    else {
      findmaximalscores_withoverlaps(chain, chaininfo, fragments,
                                     num_of_fragments, max_gap_width, seqlen1,
                                     mincoverage, chainprocessor, cpinfo,
                                     overlapinfo);
      gt_free(overlapinfo);
    }
  }
  else {
    chainingboundarycases(chain, fragments, num_of_fragments);
    if (maxscore_chains ||
        ((double) GETSTOREDLENGTH(1, 0) / (double) seqlen1) >= mincoverage) {
      chainprocessor(chain, fragments, num_of_fragments, max_gap_width, cpinfo);
    }
  }
  gt_free(chaininfo);
  gt_chain_delete(chain);
}

void gt_globalchaining_max(GtFragment *fragments,
                           unsigned long num_of_fragments,
                           unsigned long max_gap_width,
                           GtChainProc chainprocessor, void *cpinfo)
{
  globalchaining_generic(true, max_gap_width, fragments, num_of_fragments,
                         GT_UNDEF_ULONG, GT_UNDEF_DOUBLE, chainprocessor,
                         cpinfo);
}

void gt_globalchaining_coverage(GtFragment *fragments,
                                unsigned long num_of_fragments,
                                unsigned long max_gap_width,
                                unsigned long seqlen1,
                                double mincoverage, GtChainProc chainprocessor,
                                void *cpinfo)
{
  gt_assert(mincoverage >= 0.0 && mincoverage <= 1.0);
  globalchaining_generic(false, max_gap_width, fragments, num_of_fragments,
                         seqlen1, mincoverage, chainprocessor, cpinfo);
}
