/*
  Copyright (c) 2017 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2017 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "weighted_lis_filter.h"

typedef struct
{
  GtUword startpos[2],
          endpos[2],
          orgidx,
          diff,
          prev,
          distance;
  GtWord score;
  float weight;
} GtAlignmentLink;

struct GtWLisFilterMatches
{
  GtAlignmentLink *table;
  GtUword nextfree,
          allocated,
          qmax;
};

static GtUword gt_wlis_filter_aligned_len(const GtAlignmentLink *match)
{
  int idx;
  GtUword aligned_len = 0;

  gt_assert(match != NULL);
  for (idx = 0; idx < 2; idx++)
  {
    aligned_len += match->endpos[idx] - match->startpos[idx] + 1;
  }
  return aligned_len;
}

GtWLisFilterMatches *gt_wlis_filter_matches_new(void)
{
  GtWLisFilterMatches *allmatches = gt_calloc(1, sizeof (*allmatches));
  return allmatches;
}

void gt_wlis_filter_matches_reset(GtWLisFilterMatches *allmatches)
{
  gt_assert(allmatches != NULL);
  allmatches->nextfree = 0;
}

void gt_wlis_filter_matches_delete(GtWLisFilterMatches *allmatches)
{
  if (allmatches != NULL)
  {
    gt_free(allmatches->table);
    gt_free(allmatches);
  }
}

void gt_wlis_filter_matches_add(GtWLisFilterMatches *allmatches,
                                GtUword s_start, GtUword s_end,
                                GtUword q_start, GtUword q_end,
                                GtUword distance)
{
  gt_assert(allmatches);
  GtUword aligned_len;
  float prob_id;
  GtAlignmentLink *current_match;

  if (allmatches->nextfree >= allmatches->allocated)
  {
    allmatches->allocated = allmatches->allocated * 1.2 + 256;
    allmatches->table = gt_realloc(allmatches->table,
                             allmatches->allocated*sizeof (*allmatches->table));
  }
  current_match = allmatches->table + allmatches->nextfree;
  current_match->startpos[0] = s_start;
  current_match->startpos[1] = q_start;
  current_match->endpos[0] = s_end;
  current_match->endpos[1] = q_end;
  current_match->distance = distance;
  aligned_len = gt_wlis_filter_aligned_len(current_match);
  prob_id = (float) (aligned_len - 2 * distance)/aligned_len;
  current_match->weight = 10000.0 * prob_id * prob_id;
  current_match->orgidx = allmatches->nextfree++;
  if (s_end > allmatches->qmax)
  {
    allmatches->qmax = s_end;
  }
}

static int gt_alignment_link_compare(const void *vlinka, const void *vlinkb)
{
  const GtAlignmentLink *linka = (const GtAlignmentLink *) vlinka;
  const GtAlignmentLink *linkb = (const GtAlignmentLink *) vlinkb;

  if (linka->startpos[1] < linkb->startpos[1])
  {
    return -1;
  }
  if (linka->startpos[1] > linkb->startpos[1])
  {
    return 1;
  }

  if (((linka->endpos[1] - linka->startpos[1]) * linka->weight) >
      ((linkb->endpos[1] - linkb->startpos[1]) * linkb->weight))
  {
    return -1;
  }
  return 1;
}

#define GT_WLIS_FILTER_UNDEF(ABM) (ABM)->nextfree

/* leave only the alignments which form the longest mutually consistent set */
static GtUword gt_filter_apply(GtWLisFilterMatches *allmatches)
{
  GtUword bestchain_end = 0,
          maxscore = 0,
          mindiff = GT_UWORD_MAX,
          rightmatch,
          len0, len1, len;

  gt_assert(allmatches && allmatches->nextfree);
  len0 = allmatches->table[0].endpos[0] - allmatches->table[0].startpos[0]+1;
  len1 = allmatches->table[0].endpos[1] - allmatches->table[0].startpos[1]+1;
  len = MIN(len0, len1);

  allmatches->table[0].score = allmatches->table[0].weight * len;
  allmatches->table[0].diff = 0;
  allmatches->table[0].prev = GT_WLIS_FILTER_UNDEF(allmatches);

  for (rightmatch=1UL; rightmatch < allmatches->nextfree; rightmatch++)
  {
    GtUword leftmatch;

    len0 = allmatches->table[rightmatch].endpos[0] -
           allmatches->table[rightmatch].startpos[0] + 1;
    len1 = allmatches->table[rightmatch].endpos[1] -
           allmatches->table[rightmatch].startpos[1] + 1;
    len = MIN(len0, len1);

    allmatches->table[rightmatch].score
      = allmatches->table[rightmatch].weight * len;
    allmatches->table[rightmatch].diff = 0;
    allmatches->table[rightmatch].prev = GT_WLIS_FILTER_UNDEF(allmatches);

    for (leftmatch=0; leftmatch<rightmatch; leftmatch++)
    {
      GtUword diff, seq;
      GtWord overlap, score,
              ov0 = 0,
              ov1 = 0;

      /* handle overlaps */
      if (allmatches->table[leftmatch].endpos[0] >=
          allmatches->table[rightmatch].startpos[0])
      {
        ov0 = allmatches->table[leftmatch].endpos[0] -
              allmatches->table[rightmatch].startpos[0]+1;
      }

      if (allmatches->table[leftmatch].endpos[1] >=
          allmatches->table[rightmatch].startpos[1])
      {
        ov1 = allmatches->table[leftmatch].endpos[1] -
              allmatches->table[rightmatch].startpos[1]+1;
      }
      overlap = MAX(ov0, ov1);

      /* calculate score */
      score = allmatches->table[leftmatch].score +
              ((GtWord) len - overlap) * allmatches->table[rightmatch].weight;

      /* calculate Alignment gap */
      diff = allmatches->table[leftmatch].diff;
      for (seq = 0; seq < 2; seq++)
      {
        if (allmatches->table[leftmatch].startpos[seq] <
            allmatches->table[rightmatch].startpos[seq])
        {
          diff += labs(allmatches->table[leftmatch].endpos[seq] -
                       allmatches->table[rightmatch].startpos[seq]);
        }
        else
        {
          diff += labs(allmatches->table[rightmatch].endpos[seq] -
                       allmatches->table[leftmatch].startpos[seq]);
        }
      }

      /* update score */
      if (score > allmatches->table[rightmatch].score ||
         (score == allmatches->table[rightmatch].score &&
          diff < allmatches->table[rightmatch].diff))
      {
        allmatches->table[rightmatch].score = score;
        allmatches->table[rightmatch].diff = diff;
        allmatches->table[rightmatch].prev = leftmatch;

        /* update overall best score */
        if (score > maxscore || (score == maxscore && diff < mindiff))
        {
          maxscore = score;
          mindiff = diff;
          bestchain_end = rightmatch;
        }
      }
    }
  }
  return bestchain_end;
}

void gt_wlis_filter_evaluate(GtArrayGtUword *chain,
                             GtUword *sum_distance_chain,
                             GtUword *sum_aligned_len_chain,
                             GtWLisFilterMatches *allmatches,
                             bool forward)
{
  GtUword bestchain_idx, idx;
  GtUword *fwd, *bck;

  gt_assert((chain == NULL && sum_distance_chain != NULL
                           && sum_aligned_len_chain != NULL) ||
            (chain != NULL && sum_distance_chain == NULL
                           && sum_aligned_len_chain == NULL));
  /* adjust coordinates */
  if (!forward)
  {
    for (idx = 0; idx < allmatches->nextfree; idx++)
    {
      allmatches->table[idx].startpos[1]
        = allmatches->qmax - allmatches->table[idx].startpos[1];
      allmatches->table[idx].endpos[1]
        = allmatches->qmax - allmatches->table[idx].endpos[1];
    }
  }
  /* sort by query seuqence */
  qsort(allmatches->table,(size_t) allmatches->nextfree,
        sizeof *allmatches->table, gt_alignment_link_compare);

  /* call filter algorithm */
  bestchain_idx = gt_filter_apply(allmatches);

  /* get the chain by backtracing */
  do {
    if (chain != NULL)
    {
      GT_STOREINARRAY(chain,GtUword,chain->allocatedGtUword * 0.2 + 256,
                      allmatches->table[bestchain_idx].orgidx);
    } else
    {
      gt_assert(sum_distance_chain != NULL && sum_aligned_len_chain != NULL);
      *sum_distance_chain += allmatches->table[bestchain_idx].distance;
      *sum_aligned_len_chain +=
         gt_wlis_filter_aligned_len(allmatches->table + bestchain_idx);
    }
    bestchain_idx = allmatches->table[bestchain_idx].prev;
  } while (bestchain_idx != GT_WLIS_FILTER_UNDEF(allmatches));

  if (chain != NULL)
  {
  /* invert the order */
    for (fwd = chain->spaceGtUword, bck = fwd + chain->nextfreeGtUword - 1;
         fwd < bck; fwd++, bck--)
    {
      GtUword tmp = *fwd;
      *fwd = *bck;
      *bck = tmp;
    }
  }
}
