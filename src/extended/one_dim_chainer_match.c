/*
  Copyright (c) 2015 Fabian Sobanski <0sobansk@informatik.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include "extended/one_dim_chainer_match.h"
#include "extended/priority_queue.h"
#include "match/seed-extend-iter.h"

void gt_1d_chainer_match_delete(GtOneDimChainerMatch *match)
{
  GtOneDimChainerMatch *tmp = match->suc;
  gt_free(match);
  match = tmp;
  while (match != NULL)
  {
    --match->refcount;
    tmp = match->suc;
    if (match->refcount == 0)
    {
      gt_free(match);
    } else
    {
      break;
    }
    match = tmp;
  }
}

/* Starts at a given <match>, reduces its reference counter, and if there are no
references left the <match> is deleted from memory. */
void gt_1d_chainer_decr_refcount(GtOneDimChainerMatch *match)
{
  if (match != NULL)
  {
    --match->refcount;
    if (match->refcount == 0)
    {
      gt_1d_chainer_match_delete(match);
    }
  }
}

/* Increases the reference counter for a given <match>. */
void gt_1d_chainer_incr_refcount(GtOneDimChainerMatch *match)
{
  if (match != NULL)
  {
    ++match->refcount;
  }
}

GtOneDimChainerMatch* gt_1d_chainer_match_new(
    GtQuerymatch *querymatchptr, GtOneDimChainerMatch *maxchainstart,
    GtUword chainweight)
{
  GtOneDimChainerMatch *match = gt_malloc(sizeof *match);
  match->refcount = 1;
  match->querymatch = querymatchptr;
  match->suc = maxchainstart;
  gt_1d_chainer_incr_refcount(maxchainstart);
  match->chainweight = chainweight;

  return match;
}

/* Calculates and returns the weight for a match from its <querystart> and
<queryend> position. */
GtUword gt_1d_chainer_get_weight(GtUword querystart, GtUword queryend,
    GtUword dist)
{
  return gt_querymatch_distance2score(dist, queryend - querystart);
}

/* Compare function for the querystarts of 2 given matches <match1> and
<match2>. Returns 1 if the end of <match1> is smaller (i.e. its priority is
greater), 0 if they are equal, and -1 if the end of <match2> is smaller (i.e.
the priority of <match2> is greater). This method is used for creating the
priority queue. */
int gt_1d_chainer_compare_starts(const void *match1, const void *match2)
{
  GtUword start1 = gt_querymatch_querystart(
      ((GtOneDimChainerMatch*) match1)->querymatch);
  GtUword start2 = gt_querymatch_querystart(
      ((GtOneDimChainerMatch*) match2)->querymatch);

  if (start1 > start2) {
    return -1;
  } else if (start1 == start2) {
    return 0;
  } else {
    return 1;
  }
}

int gt_1d_chainer_calc_chain(GtUword overlap, GtOneDimChainerMatch *chainstart,
                             GtSeedextendMatchIterator *semi,
                             GT_UNUSED GtError *err)
{
  int had_err = 0;

  /* Set up the priority queue for match iteration */
  GtUword maxnumofelements = 15;
  GtPriorityQueue *pq = gt_priority_queue_new(gt_1d_chainer_compare_starts,
         maxnumofelements);
  if (pq == NULL)
  {
    return -1;
  }

  /* Initialize important chaining variables */
  GtUword chainweight = 0;
  GtOneDimChainerMatch *match = NULL;
  GtOneDimChainerMatch *maxchainstart = NULL;
  uint64_t lastseqnum = 0;

  /* Begin iteration over given matches */
  while (true)
  {
    /* Get new match from iterator */
    GtQuerymatch *querymatchptr = gt_seedextend_match_iterator_next(semi);

    /* Iterate over priority queue */
    while (!gt_priority_queue_is_empty(pq))
    {
      GtOneDimChainerMatch *candidatematch =
        (GtOneDimChainerMatch*) gt_priority_queue_find_min(pq);
      GtQuerymatch *candidatesem = candidatematch->querymatch;
      GtUword candidateend = gt_querymatch_querystart(candidatesem) +
                             gt_querymatch_querylen(candidatesem);
      GtOneDimChainerMatch *candidatesuc = candidatematch->suc;
      if (candidatesuc != NULL &&
          lastseqnum == gt_querymatch_queryseqnum(candidatesuc->querymatch))
      {
        candidateend = fmin(candidateend,
                            gt_querymatch_querystart(candidatesuc->querymatch));
      }

      /* Continue only if <querymatchptr> is defined, we are comparing
         matches from the same genome, and the candidate does not overlap
         more than allowed with the match from the iterator. */
      if (querymatchptr != NULL &&
          lastseqnum == gt_querymatch_queryseqnum(querymatchptr) &&
          fmin(gt_querymatch_querystart(candidatesem) + overlap,
               candidateend - 1) < gt_querymatch_querystart(querymatchptr) +
                                   gt_querymatch_querylen(querymatchptr))
      {
        break;
      } else
      {
        /* Remove candidate from priority queue */
        gt_priority_queue_extract_min(pq);

        /* If the weight of the candidate added to the chain increases
           the maximum chain weight, extend the chain. */
        GtUword candidatechainweight = candidatematch->chainweight +
          gt_1d_chainer_get_weight(gt_querymatch_querystart(candidatesem),
              candidateend, gt_querymatch_distance(candidatesem));
        if (chainweight < candidatechainweight)
        {
          chainweight = candidatechainweight;
          gt_1d_chainer_decr_refcount(maxchainstart);
          maxchainstart = candidatematch;
        } else
        {
          gt_1d_chainer_decr_refcount(candidatematch);
        }
      }
    }

    if (querymatchptr == NULL)
    {
      break;
    }

    /* Update the sequence ID if necessary */
    if (lastseqnum != gt_querymatch_queryseqnum(querymatchptr))
    {
      lastseqnum = gt_querymatch_queryseqnum(querymatchptr);
    }

    /* Create a new chain match object from the given iterator match */
    gt_1d_chainer_decr_refcount(match);
    match = gt_1d_chainer_match_new(querymatchptr, maxchainstart, chainweight);
    if (match == NULL)
    {
      return -1; /* We are not able to free previously allocated matches. :( */
    }

    /* Extend priority queue in case capacity is not enough.
       Almost never happens.*/
    if (gt_priority_queue_is_full(pq))
    {
      maxnumofelements *= 2;
      GtPriorityQueue *newpq = gt_priority_queue_new(
          gt_1d_chainer_compare_starts, maxnumofelements);
      if (newpq == NULL)
      {
        return -1; /* We are not able to free previously allocated matches. */
      }
      while (!gt_priority_queue_is_empty(pq))
      {
        gt_priority_queue_add(newpq, gt_priority_queue_extract_min(pq));
      }
      gt_priority_queue_delete(pq);
      pq = newpq;
    }
    gt_1d_chainer_incr_refcount(match);
    gt_priority_queue_add(pq, match);
  }

  *chainstart = *maxchainstart;
  gt_1d_chainer_incr_refcount(chainstart->suc);

  /* Delete priority queue after chain has been created */
  gt_priority_queue_delete(pq);
  gt_1d_chainer_decr_refcount(match);
  gt_1d_chainer_decr_refcount(maxchainstart);

  return had_err;
}
