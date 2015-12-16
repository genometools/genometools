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
  GtOneDimChainerMatch *tmp = match->prec;
  gt_free(match);
  match = tmp;
  while (match != NULL)
  {
    --match->refcount;
    tmp = match->prec;
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
    GtQuerymatch *querymatchptr, GtOneDimChainerMatch *maxchainend,
    GtUword chainweight)
{
  GtOneDimChainerMatch *match = gt_malloc(sizeof *match);
  match->refcount = 1;
  match->seqnum = gt_querymatch_queryseqnum(querymatchptr);
  match->start = gt_querymatch_querystart(querymatchptr);
  match->end = match->start + gt_querymatch_querylen(querymatchptr);
  match->prec = maxchainend;
  gt_1d_chainer_incr_refcount(maxchainend);
  match->chainweight = chainweight;
  match->dist = gt_querymatch_distance(querymatchptr);

  return match;
}

/* Calculates and returns the weight for a match from its <start> and <end>
position.
The weight equals the length of the match. */
GtUword gt_1d_chainer_get_weight(GtUword start, GtUword end,
    GtUword dist)
{
  return gt_querymatch_distance2score(dist, end - start);
}

/* Compare function for the ends of 2 given matches <match1> and <match2>.
Returns -1 if <match1> is smaller, 0 if they are equal, and 1 if <match2> end
is smaller. This method is used for creating the priority queue. */
int gt_1d_chainer_compare_ends(const void *match1, const void *match2)
{
  const GtOneDimChainerMatch *m1 = (GtOneDimChainerMatch*) match1;
  const GtOneDimChainerMatch *m2 = (GtOneDimChainerMatch*) match2;
  if (m1->end < m2->end) {
    return -1;

  } else if (m1->end == m2->end) {
    return 0;

  } else {
    return 1;
  }
}

int gt_1d_chainer_calc_chain(const GtStr *matchfilename, GtUword overlap,
                             GtOneDimChainerMatch *chainend, GtError *err)
{
  int had_err = 0;

  /* Set up match iterator */
  GtSeedextendMatchIterator *semi
       = gt_seedextend_match_iterator_new(matchfilename, err);

  if (semi == NULL)
  {
    return -1;
  }

  /* Set up the priority queue for match iteration */
  GtUword maxnumofelements = 15;
  GtPriorityQueue *pq = gt_priority_queue_new(gt_1d_chainer_compare_ends,
         maxnumofelements);
  if (pq == NULL)
  {
    return -1;
  }

  /* Initialize important chaining variables */
  GtUword chainweight = 0;
  GtOneDimChainerMatch *match = NULL;
  GtOneDimChainerMatch *maxchainend = NULL;
  uint64_t lastseqnum = 0;

  /* Begin iteration over given matches */
  while (true)
  {
    /* Get new match from iterator */
    GtQuerymatch *querymatchptr = gt_seedextend_match_iterator_next(semi);

    /* Iterate over priority queue */
    while (!gt_priority_queue_is_empty(pq))
    {
      /* Match with the minimum start position is a candidate */
      GtOneDimChainerMatch *candidatematch =
        (GtOneDimChainerMatch*) gt_priority_queue_find_min(pq);
      GtUword start = candidatematch->start;

      if (maxchainend != NULL)
      {
        start = fmax(start, maxchainend->end);
      }

      /* Continue only if <querymatchptr> is defined, we are comparing
         matches from the same genome, and the candidate does not overlap
         more than allowed with the match from the iterator. */
      if (querymatchptr != NULL &&
          lastseqnum == gt_querymatch_queryseqnum(querymatchptr) &&
          fmax(candidatematch->end - overlap, start + 1) >
          gt_querymatch_querystart(querymatchptr))
      {
        break;
      } else
      {
        /* Remove candidate from priority queue */
        gt_priority_queue_extract_min(pq);

        /* If the weight of the candidate added to the chain increases
           the maximum chain weight, extend the chain. */
        if (chainweight < candidatematch->chainweight +
            gt_1d_chainer_get_weight(start, candidatematch->end,
              candidatematch->dist))
        {
          chainweight = candidatematch->chainweight +
            gt_1d_chainer_get_weight(start, candidatematch->end,
                candidatematch->dist);
          gt_1d_chainer_decr_refcount(maxchainend);
          maxchainend = candidatematch;
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
    match = gt_1d_chainer_match_new(querymatchptr, maxchainend, chainweight);
    if (match == NULL)
    {
      return -1; /* We are not able to free previously allocated matches. :( */
    }

    /* Extend priority queue in case capacity is not enough.
       Almost never happens.*/
    if (gt_priority_queue_is_full(pq))
    {
      maxnumofelements *= 2;
      GtPriorityQueue *newpq = gt_priority_queue_new(gt_1d_chainer_compare_ends,
          maxnumofelements);
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

  *chainend = *maxchainend;
  gt_1d_chainer_incr_refcount(chainend->prec);

  /* Delete priority queue and iterator after chain has been created */
  gt_priority_queue_delete(pq);
  gt_seedextend_match_iterator_delete(semi);
  gt_1d_chainer_decr_refcount(match);
  gt_1d_chainer_decr_refcount(maxchainend);

  return had_err;
}
