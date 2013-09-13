/*
  Copyright (c) 2007, 2009 Gordon Gremme <gordon@gremme.org>
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
  GtWord maxscore;
  GtUword maxfragnum;
  bool defined;
} Maxfragvalue;

typedef struct {
  GtUword previousinchain;  /* previous index in chain */
  GtWord score;             /* score of highest scoring chain ending here,
                               computed */
} GtChaininfo;

/*
  The following structure is used to store additional information for every
  fragment when chains containing overlaps are computed.
*/

typedef struct {
  bool active;                     /* fragment is active, i.e., it is an end
                                      fragment */
  GtUword startofchain,      /* fragment number starting this chain */
                dim1lengthofchain, /* length of chain in the first dimension
                                      (necessary to compute the coverage) */
                chainarray;        /* contains computed chains:
                                      [startofchain] -> end fragment number */
} Overlapinfo;

static void initoverlapinfo(Overlapinfo *overlapinfo,
                            GtFragment *fragments,
                            GtUword num_of_fragments)
{
  GtUword i;

  for (i = 0; i < num_of_fragments; i++) {
    overlapinfo[i].active            = true;
    overlapinfo[i].startofchain      = i;
    overlapinfo[i].dim1lengthofchain = GETSTOREDLENGTH(1, i);
    overlapinfo[i].chainarray        = UNDEFPREVIOUS;
  }
}

static bool colinearfragments(GtFragment *fragments,
                              GtUword i, GtUword j)
{
  if (GETSTOREDSTARTPOINT(1, i) < GETSTOREDSTARTPOINT(1, j) &&
      GETSTOREDENDPOINT(1, i)   < GETSTOREDENDPOINT(1, j)   &&
      GETSTOREDSTARTPOINT(2, i) < GETSTOREDSTARTPOINT(2, j) &&
      GETSTOREDENDPOINT(2, i)   < GETSTOREDENDPOINT(2, j)) {
    return true;
  }
  return false;
}

static GtWord overlapcost(GtFragment *fragments,
                        GtUword i, GtUword j)
{
  GtUword overlaplength = 0;

  /* add overlap in first dimension */
  if (GETSTOREDSTARTPOINT(1, j) <= GETSTOREDENDPOINT(1, i))
    overlaplength += GETSTOREDENDPOINT(1, i) - GETSTOREDSTARTPOINT(1, j) + 1;

  /* add overlap in second dimension */
  if (GETSTOREDSTARTPOINT(2, j) <= GETSTOREDENDPOINT(2, i))
    overlaplength += GETSTOREDENDPOINT(2, i) - GETSTOREDSTARTPOINT(2, j) + 1;

  gt_log_log("overlap total  (#"GT_WU", #"GT_WU")="GT_WU"", i, j,
             overlaplength);

  return (GtWord) overlaplength;
}

static void chainingboundarycases(GtChain *chain,
                                  GtFragment *fragments,
                                  GtUword num_of_fragments)
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
                                   GtUword num_of_fragments,
                                   GtUword retracestart)
{
  GtUword fragnum, idx, lengthofchain;

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
                                GtUword max_gap_width,
                                GtUword leftfrag,
                                GtUword rightfrag)
{
  GtUword gapwidth, startpoint, endpoint;

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
                                     GtUword max_gap_width,
                                     GtFragment *fragments,
                                     GtUword num_of_fragments,
                                     Overlapinfo *overlapinfo)
{
  GtUword previous, leftfrag, rightfrag, overlaplength;
  GtWord weightright, score;
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
                                     GtUword num_of_fragments,
                                     GtUword currentfrag)
{
  if (currentfrag == num_of_fragments - 1)
    return true;
  if (chaininfo[currentfrag+1].previousinchain != currentfrag)
    return true;
  if (chaininfo[currentfrag+1].score < chaininfo[currentfrag].score)
    return true;
  return false;
}

static bool retrievemaximalscore(GtWord *maxscore, GtChaininfo *chaininfo,
                                 GtUword num_of_fragments)
{
  GtUword i;
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
                                   GtUword num_of_fragments,
                                   GtUword max_gap_width,
                                   GtChain *chain, GtWord minscore,
                                   GtChainProc chainprocessor,
                                   void *cpinfo)
{
  GtUword i;

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
                              GtUword num_of_fragments,
                              GtUword max_gap_width,
                              GtChainProc chainprocessor, void *cpinfo)
{
  GtWord minscore;
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
                                           GtUword num_of_fragments,
                                           GtUword max_gap_width,
                                           GtUword seqlen1,
                                           double mincoverage,
                                           GtChainProc chainprocessor,
                                           void *cpinfo,
                                           Overlapinfo *overlapinfo)
{
  GtUword i, startfrag;
  GtArray *startfragments;

  gt_assert(seqlen1 != GT_UNDEF_UWORD);
  gt_assert(mincoverage != GT_UNDEF_DOUBLE);
  startfragments = gt_array_new(sizeof (GtUword));

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
    startfrag = *(GtUword*) gt_array_get(startfragments, i);
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

static void log_fragments(GtFragment *fragments, GtUword num_of_fragments)
{
  GtUword i;
  gt_log_log("show chaining fragments");
  for (i = 0; i < num_of_fragments; i++) {
    GtFragment *frag = fragments + i;
    gt_log_log("#"GT_WU": s1="GT_WU", s1="GT_WU", l1="GT_WU", s2="GT_WU", "
               "e2="GT_WU", l2="GT_WU", w="GT_WU, i, frag->startpos1,
               frag->endpos1, frag->endpos1 - frag->startpos1 + 1,
               frag->startpos2, frag->endpos2,
               frag->endpos2 - frag->startpos2 + 1, frag->weight);
  }
}

static void globalchaining_generic(bool maxscore_chains,
                                   GtUword max_gap_width,
                                   GtFragment *fragments,
                                   GtUword num_of_fragments,
                                   GtUword seqlen1, double mincoverage,
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
                           GtUword num_of_fragments,
                           GtUword max_gap_width,
                           GtChainProc chainprocessor, void *cpinfo)
{
  globalchaining_generic(true, max_gap_width, fragments, num_of_fragments,
                         GT_UNDEF_UWORD, GT_UNDEF_DOUBLE, chainprocessor,
                         cpinfo);
}

void gt_globalchaining_coverage(GtFragment *fragments,
                                GtUword num_of_fragments,
                                GtUword max_gap_width,
                                GtUword seqlen1,
                                double mincoverage, GtChainProc chainprocessor,
                                void *cpinfo)
{
  gt_assert(mincoverage >= 0.0 && mincoverage <= 1.0);
  globalchaining_generic(false, max_gap_width, fragments, num_of_fragments,
                         seqlen1, mincoverage, chainprocessor, cpinfo);
}
