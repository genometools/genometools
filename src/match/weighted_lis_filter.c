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

#define GT_UNDEFPREVIOUS allmatches->nextfree
#define ADD_SIZE 10

#include "core/assert_api.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "weighted_lis_filter.h"

GtFilter *gt_filter_new()
{
  GtFilter *filter = gt_calloc(1, sizeof (*filter));
  filter->allocated = ADD_SIZE;
  filter->chain = gt_malloc(sizeof (*filter->chain) * filter->allocated);

  return filter;
}

void gt_filter_reset(GtFilter *filter)
{
  gt_assert(filter);
  filter->nextfree = 0;
}

void gt_filter_delete(GtFilter *filter)
{
  if (filter != NULL)
  {
    gt_free(filter->chain);
    gt_free(filter);
  }
}

typedef struct
{
  GtUword startpos[2],
          endpos[2],
          orgidx,
          diff,
          prev;
  GtWord  score;
  float   weight;

}GtAlignmentLink;

struct GtAllMatches
{
  GtAlignmentLink *table;
  GtUword nextfree,
          allocated,
          qmax;
  bool reverse;
};

GtAllMatches *gt_filter_all_matches_new(bool reverse)
{
  GtAllMatches *allmatches;
  allmatches = gt_calloc(1, sizeof (*allmatches));
  allmatches->reverse = reverse;

  return allmatches;
}

void gt_filter_all_matches_reset(GtAllMatches *allmatches, bool reverse)
{
  gt_assert(allmatches);
  allmatches->nextfree = 0;
  allmatches->reverse = reverse;
}

void gt_filter_all_matches_delete(GtAllMatches *allmatches)
{
  if (allmatches != NULL)
  {
    gt_free(allmatches->table);
    gt_free(allmatches);
  }
}

void gt_filter_all_matches_add(GtAllMatches *allmatches,
                               GtUword s_start, GtUword s_end,
                               GtUword q_start, GtUword q_end, float weight)
{
  gt_assert(allmatches);

  if (allmatches->nextfree >= allmatches->allocated)
  {
    allmatches->allocated += ADD_SIZE;
    allmatches->table = gt_realloc(allmatches->table,
                             allmatches->allocated*sizeof (*allmatches->table));
  }

  allmatches->table[allmatches->nextfree].startpos[0] = s_start;
  allmatches->table[allmatches->nextfree].startpos[1] = q_start;
  allmatches->table[allmatches->nextfree].endpos[0] = s_end;
  allmatches->table[allmatches->nextfree].endpos[1] = q_end;
  allmatches->table[allmatches->nextfree].weight = weight;
  allmatches->table[allmatches->nextfree].orgidx = allmatches->nextfree;
  allmatches->nextfree++;

  if (s_end > allmatches->qmax)
    allmatches->qmax = s_end;

}

static int querycomparefunc(const void *linka, const void *linkb)
{
  if (((const GtAlignmentLink *) linka)->startpos[1] <
      ((const GtAlignmentLink *) linkb)->startpos[1])
  {
    return -1;
  }
  if (((const GtAlignmentLink *) linka)->startpos[1] >
      ((const GtAlignmentLink *) linkb)->startpos[1])
  {
    return 1;
  }

   if (((((const GtAlignmentLink *) linka)->endpos[1]-
         ((const GtAlignmentLink *) linka)->startpos[1]) *
         ((const GtAlignmentLink *) linka)->weight ) >
         ((((const GtAlignmentLink *) linkb)->endpos[1]-
         ((const GtAlignmentLink *) linkb)->startpos[1]) *
         ((const GtAlignmentLink *) linkb)->weight ))
     return -1;
   else
     return 1;

  return 0;
}

/* leave only the alignments which form the longest mutually consistent set */
GtUword gt_filter_apply(GtAllMatches *allmatches)
{
  gt_assert(allmatches && allmatches->nextfree);

  GtUword bestchain_end = 0,
          maxscore = 0,
          mindiff = GT_UWORD_MAX,
          rightmatch,
          len0, len1, len;

  len0 = allmatches->table[0].endpos[0] - allmatches->table[0].startpos[0]+1;
  len1 = allmatches->table[0].endpos[1] - allmatches->table[0].startpos[1]+1;
  len = MIN(len0, len1);

  /* initialize */
  allmatches->table[0].score = allmatches->table[0].weight * len;
  allmatches->table[0].diff = 0;
  allmatches->table[0].prev = GT_UNDEFPREVIOUS;

  for (rightmatch=1UL; rightmatch < allmatches->nextfree; rightmatch++)
  {

    GtUword leftmatch;

    len0 = allmatches->table[rightmatch].endpos[0] -
           allmatches->table[rightmatch].startpos[0] + 1;
    len1 = allmatches->table[rightmatch].endpos[1] -
           allmatches->table[rightmatch].startpos[1] + 1;
    len = MIN(len0, len1);

    allmatches->table[rightmatch].score =
                                   allmatches->table[rightmatch].weight * len;
    allmatches->table[rightmatch].diff = 0;
    allmatches->table[rightmatch].prev = GT_UNDEFPREVIOUS;

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
                  ((GtWord)len-overlap) * allmatches->table[rightmatch].weight;

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

void gt_filter_call(GtAllMatches *allmatches, GtFilter *filter)
{
  GtUword bestchain_idx, idx, tmp;

  /* adjust coordinates */
  if (allmatches->reverse)
  {
    for (idx = 0; idx < allmatches->nextfree; idx++)
    {
      allmatches->table[idx].startpos[1] =
                      allmatches->qmax -allmatches->table[idx].startpos[1];
      allmatches->table[idx].endpos[1] =
                      allmatches->qmax -allmatches->table[idx].endpos[1];
    }
  }
  /* sort by query seuqence */
  qsort(allmatches->table,(size_t) allmatches->nextfree,
        sizeof (*allmatches->table), querycomparefunc);

  /* call filter algorithm */
  bestchain_idx = gt_filter_apply(allmatches);

  /* get the chain by backtracing */
  do {

    if (filter->nextfree >= filter->allocated)
    {
      filter->allocated += ADD_SIZE;
      filter->chain = gt_realloc(filter->chain,
                                 filter->allocated * sizeof (*filter->chain));
    }
    filter->chain[filter->nextfree] = allmatches->table[bestchain_idx].orgidx;
    filter->nextfree += 1;
    bestchain_idx = allmatches->table[bestchain_idx].prev;

  } while (bestchain_idx != GT_UNDEFPREVIOUS);

  /* invert the order */
  for (idx = 0; idx < filter->nextfree/2; idx++)
  {
    tmp = filter->chain[idx];
    filter->chain[idx] = filter->chain[filter->nextfree-idx-1];
    filter->chain[filter->nextfree-1-idx] = tmp;
  }
}
