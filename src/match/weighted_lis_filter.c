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
  GtWord score;
  GtUword startpos[2],
          endpos[2],
          diff,
          prev;
  union
  {
    GtUword original_index, /* Only needed when matches are output */
            distance; /* not needed when matches are output */
  } oi_di;
  float weight;
} GtWlisItem;

static GtUword gt_wlis_filter_aligned_len(const GtWlisItem *match)
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

GT_DECLAREARRAYSTRUCT(GtWlisItem);

struct GtWLisFilterMatches
{
  GtArrayGtWlisItem items;
  GtUword qmax;
};

GtWLisFilterMatches *gt_wlis_filter_matches_new(void)
{
  GtWLisFilterMatches *wlismatches = gt_malloc(sizeof *wlismatches);

  GT_INITARRAY(&wlismatches->items,GtWlisItem);
  wlismatches->qmax = 0;
  return wlismatches;
}

void gt_wlis_filter_matches_reset(GtWLisFilterMatches *wlismatches)
{
  gt_assert(wlismatches != NULL);
  wlismatches->items.nextfreeGtWlisItem = 0;
  wlismatches->qmax = 0;
}

void gt_wlis_filter_matches_delete(GtWLisFilterMatches *wlismatches)
{
  if (wlismatches != NULL)
  {
    GT_FREEARRAY(&wlismatches->items,GtWlisItem);
    gt_free(wlismatches);
  }
}

void gt_wlis_filter_matches_add(GtWLisFilterMatches *wlismatches,
                                GtUword s_start, GtUword s_end,
                                GtUword q_start, GtUword q_end,
                                GtUword distance,
                                bool store_querymatch)
{
  GtUword aligned_len;
  float prob_id;
  GtWlisItem *current_match;

  gt_assert(wlismatches != NULL);
  GT_GETNEXTFREEINARRAY(current_match,&wlismatches->items,GtWlisItem,
                        wlismatches->items.
                                     allocatedGtWlisItem * 0.2 + 256);
  current_match->startpos[0] = s_start;
  current_match->startpos[1] = q_start;
  current_match->endpos[0] = s_end;
  current_match->endpos[1] = q_end;

  if (store_querymatch)
  {
    gt_assert(wlismatches->items.nextfreeGtWlisItem > 0);
    current_match->oi_di.original_index
      = wlismatches->items.nextfreeGtWlisItem - 1;
  } else
  {
    current_match->oi_di.distance = distance;
  }
  aligned_len = gt_wlis_filter_aligned_len(current_match);
  prob_id = (float) (aligned_len - 2 * distance)/aligned_len;
  current_match->weight = prob_id * prob_id;


  if (q_end > wlismatches->qmax)
  {
    wlismatches->qmax = q_end;
  }
}

static int gt_alignment_link_compare(const void *vlinka, const void *vlinkb)
{
  const GtWlisItem *linka = (const GtWlisItem *) vlinka;
  const GtWlisItem *linkb = (const GtWlisItem *) vlinkb;

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

#define GT_WLIS_FILTER_UNDEF(ABM) (ABM)->items.nextfreeGtWlisItem
#define GT_WLIS_ACC(IDX)          wlismatches->items.spaceGtWlisItem[IDX]

/* leave only the alignments which form the longest mutually consistent set */
static GtUword gt_filter_apply(GtWLisFilterMatches *wlismatches)
{
  GtUword bestchain_end = 0,
          maxscore = 0,
          mindiff = GT_UWORD_MAX,
          rightmatch,
          len0, len1, len;

  gt_assert(wlismatches != NULL &&
            wlismatches->items.nextfreeGtWlisItem > 0);
  len0 = GT_WLIS_ACC(0).endpos[0] - GT_WLIS_ACC(0).startpos[0]+1;
  len1 = GT_WLIS_ACC(0).endpos[1] - GT_WLIS_ACC(0).startpos[1]+1;
  len = MIN(len0, len1);

  GT_WLIS_ACC(0).score = GT_WLIS_ACC(0).weight * len;
  GT_WLIS_ACC(0).diff = 0;
  GT_WLIS_ACC(0).prev = GT_WLIS_FILTER_UNDEF(wlismatches);

  for (rightmatch=1UL; rightmatch < wlismatches->items.nextfreeGtWlisItem;
       rightmatch++)
  {
    GtUword leftmatch;

    len0 = GT_WLIS_ACC(rightmatch).endpos[0] -
           GT_WLIS_ACC(rightmatch).startpos[0] + 1;
    len1 = GT_WLIS_ACC(rightmatch).endpos[1] -
           GT_WLIS_ACC(rightmatch).startpos[1] + 1;
    len = MIN(len0, len1);

    GT_WLIS_ACC(rightmatch).score = GT_WLIS_ACC(rightmatch).weight * len;
    GT_WLIS_ACC(rightmatch).diff = 0;
    GT_WLIS_ACC(rightmatch).prev = GT_WLIS_FILTER_UNDEF(wlismatches);

    for (leftmatch=0; leftmatch<rightmatch; leftmatch++)
    {
      int dim;
      GtUword diff, overlap, ovtab[2] = {0};
      GtWord score;
      
      /* handle overlaps and calculate Alignment gap */
      diff = GT_WLIS_ACC(leftmatch).diff;
      for (dim = 0; dim < 2; dim++)
      {
        if (GT_WLIS_ACC(leftmatch).endpos[dim] >=
            GT_WLIS_ACC(rightmatch).startpos[dim])
        {
          ovtab[dim] = GT_WLIS_ACC(leftmatch).endpos[dim] -
                       GT_WLIS_ACC(rightmatch).startpos[dim] + 1;
        }
        if (GT_WLIS_ACC(leftmatch).startpos[dim] <
            GT_WLIS_ACC(rightmatch).startpos[dim])
        {
          diff += labs((GtWord) GT_WLIS_ACC(leftmatch).endpos[dim] -
                       (GtWord) GT_WLIS_ACC(rightmatch).startpos[dim]);
        } else
        {
          diff += labs((GtWord) GT_WLIS_ACC(rightmatch).endpos[dim] -
                       (GtWord) GT_WLIS_ACC(leftmatch).startpos[dim]);
        }
      }
      overlap = MAX(ovtab[0], ovtab[1]);
      
      /* calculate score */
      score = GT_WLIS_ACC(leftmatch).score +
              (GtWord)(((GtWord) len - (GtWord) overlap) *
              GT_WLIS_ACC(rightmatch).weight);

      /* update score */
      if (score > GT_WLIS_ACC(rightmatch).score ||
         (score == GT_WLIS_ACC(rightmatch).score &&
          diff < GT_WLIS_ACC(rightmatch).diff))
      {
        GT_WLIS_ACC(rightmatch).score = score;
        GT_WLIS_ACC(rightmatch).diff = diff;
        GT_WLIS_ACC(rightmatch).prev = leftmatch;

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
                             GtWLisFilterMatches *wlismatches,
                             bool forward)
{
  GtUword bestchain_idx, *fwd, *bck;
  
  if (wlismatches->items.nextfreeGtWlisItem == 0)
  {
    return;
  }
  gt_assert((chain == NULL && sum_distance_chain != NULL
                           && sum_aligned_len_chain != NULL) ||
            (chain != NULL && sum_distance_chain == NULL
                           && sum_aligned_len_chain == NULL));
  /* adjust coordinates */
  if (!forward)
  {
    GtUword idx;
    for (idx = 0; idx < wlismatches->items.nextfreeGtWlisItem; idx++)
    {
      GtUword tmp = GT_WLIS_ACC(idx).startpos[1];
      GT_WLIS_ACC(idx).startpos[1]
        = wlismatches->qmax - GT_WLIS_ACC(idx).endpos[1];
      GT_WLIS_ACC(idx).endpos[1]
        = wlismatches->qmax - tmp;
    }
  }

  /* sort by query seuqence */
  qsort(wlismatches->items.spaceGtWlisItem,
        (size_t) wlismatches->items.nextfreeGtWlisItem,
        sizeof *wlismatches->items.spaceGtWlisItem,
        gt_alignment_link_compare);

  /* call filter algorithm */
  bestchain_idx = gt_filter_apply(wlismatches);

  /* get the chain by backtracing */
  do {
    if (chain != NULL)
    {
      GT_STOREINARRAY(chain,GtUword,chain->allocatedGtUword * 0.2 + 256,
                      GT_WLIS_ACC(bestchain_idx).oi_di.original_index);
    } else
    {
      gt_assert(sum_distance_chain != NULL && sum_aligned_len_chain != NULL);
      *sum_distance_chain += GT_WLIS_ACC(bestchain_idx).oi_di.distance;
      *sum_aligned_len_chain +=
         gt_wlis_filter_aligned_len(wlismatches->items.spaceGtWlisItem +
                                    bestchain_idx);
    }
    bestchain_idx = GT_WLIS_ACC(bestchain_idx).prev;
  } while (bestchain_idx != GT_WLIS_FILTER_UNDEF(wlismatches));

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
