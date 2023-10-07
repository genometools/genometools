#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <limits.h>
#ifdef WITHDISTRIBUTION
#include "distribution.h"
#endif
#include "core/assert_api.h"
#include "core/divmodmul_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "match/ft-front-generation.h"
#include "match/ft-eoplist.h"

#define GT_FT_BACKTRACEBITS 3

typedef struct
{
  uint32_t bits:GT_FT_BACKTRACEBITS,     /* combination of FT_EOP_MISMATCH
                                                           FT_EOP_INSERTION
                                                           FT_EOP_DELETION */
           lcs:(32 - GT_FT_BACKTRACEBITS); /* longest common suffix */
} GtBackreftable;

typedef struct
{
  uint8_t eopcode;
  unsigned int lcs;
} GtBacktraceFrontpath;

typedef struct
{
  GtWord diagonal,
         scoresum;
  GtUword distance,
          globaloffset,
          trimleft,
          lcs_sum,
          pathlength;
  unsigned int row, lcs;
  uint8_t trace,
          eopcode;
} GtBacktraceFrontStackelem;

typedef struct
{
  GtBacktraceFrontStackelem *space;
  GtUword nextfree, allocated;
} GtBacktraceFrontStack;

typedef uint16_t GtFrontGenerationValue;
#define GT_FRONTGENERATION_VALUE_MAX\
        ((1UL << (sizeof(GtFrontGenerationValue) * CHAR_BIT)) - 1)

typedef struct
{
  GtFrontGenerationValue trimleft_diff,
                         valid;
} GtFrontGeneration;

struct GtFrontTrace
{
  GtBackreftable *backref_table;
  GtUword maxlcs,
          backref_nextfree,
          backref_allocated;
  /* from now on all values are used for the trimming version only */
  GtFrontGeneration *gen_table;
#ifdef WITHDISTRIBUTION
  Distribution *lcs_dist, *valid_dist, *trimleft_diff_dist, *space_per_pos_dist;
#endif
  GtBacktraceFrontpath *backtracepath;
  GtBacktraceFrontStack backtracestack;
  GtUword backtracepath_allocated,
          gen_nextfree,
          gen_allocated,
          previoustrimleft;
};

GtFrontTrace *front_trace_new(void)
{
  GtFrontTrace *front_trace = gt_malloc(sizeof *front_trace);

  gt_assert(front_trace != NULL);
  /*printf("sizeof *backref_table=" GT_WU " bytes\n",
          sizeof *front_trace->backref_table);*/
  front_trace->backref_table = NULL;
  front_trace->backref_nextfree = 0;
  front_trace->backref_allocated = 0;
  front_trace->gen_table = NULL;
  front_trace->maxlcs = (1 << (32 - GT_FT_BACKTRACEBITS)) - 1;
  front_trace->gen_allocated = 0;
  front_trace->gen_nextfree = 0;
  front_trace->previoustrimleft = 0;
#ifdef WITHDISTRIBUTION
  front_trace->lcs_dist = distribution_new();
  front_trace->valid_dist = distribution_new();
  front_trace->trimleft_diff_dist = distribution_new();
  front_trace->space_per_pos_dist = distribution_new();
#endif
  front_trace->backtracepath = NULL;
  front_trace->backtracepath_allocated = 0;
  front_trace->backtracestack.nextfree = front_trace->backtracestack.allocated
                                       = 0;
  front_trace->backtracestack.space = NULL;
  return front_trace;
}

#ifdef WITHDISTRIBUTION
static size_t front_trace_size(const GtFrontTrace *front_trace)
{
  gt_assert(front_trace != NULL);
  return sizeof *front_trace->gen_table * front_trace->gen_nextfree +
         sizeof *front_trace->backref_table * front_trace->backref_nextfree;
}
#endif

void front_trace_reset(GtFrontTrace *front_trace,
                       GT_UNUSED GtUword sumseqlen)
{
#ifdef WITHDISTRIBUTION
  double space_per_pos = (double) front_trace_size(front_trace)/sumseqlen;
#endif

  gt_assert(front_trace != NULL);
#ifdef WITHDISTRIBUTION
  distribution_add(front_trace->space_per_pos_dist,
                   (GtUword) (space_per_pos * 100.0));
#endif
  front_trace->backref_nextfree = 0;
  front_trace->gen_nextfree = 0;
}

void front_trace_delete(GtFrontTrace *front_trace)
{
  if (front_trace != NULL)
  {
    /*
    distribution_show("lcs",1,front_trace->lcs_dist);
    distribution_show("valid",1,front_trace->valid_dist);
    distribution_show("trimleft_diff",1,front_trace->trimleft_diff_dist);
    distribution_show("space_per_pos",100,front_trace->space_per_pos_dist);
    */
#ifdef WITHDISTRIBUTION
    distribution_delete(front_trace->lcs_dist);
    distribution_delete(front_trace->valid_dist);
    distribution_delete(front_trace->trimleft_diff_dist);
    distribution_delete(front_trace->space_per_pos_dist);
#endif
    gt_free(front_trace->backref_table);
    gt_free(front_trace->gen_table);
    gt_free(front_trace->backtracepath);
    gt_free(front_trace->backtracestack.space);
    gt_free(front_trace);
  }
}

void front_trace_add_gen(GtFrontTrace *front_trace,GtUword trimleft,
                         GtUword valid)
{
  gt_assert (front_trace != NULL);
  if (front_trace->gen_nextfree >= front_trace->gen_allocated)
  {
    front_trace->gen_allocated
      = front_trace->gen_allocated * 1.2 + 128UL;
    front_trace->gen_table
      = gt_realloc(front_trace->gen_table,
                   sizeof *front_trace->gen_table * front_trace->gen_allocated);
    gt_assert(front_trace->gen_table != NULL);
  }
#ifdef WITHDISTRIBUTION
  distribution_add(front_trace->valid_dist,valid);
#endif
  if (front_trace->gen_nextfree > 0)
  {
    GtUword trimleft_diff;

    gt_assert(front_trace->previoustrimleft <= trimleft);
    trimleft_diff = trimleft - front_trace->previoustrimleft;
#ifdef WITHDISTRIBUTION
    distribution_add(front_trace->trimleft_diff_dist,trimleft_diff);
#endif
    gt_assert(trimleft_diff <= GT_FRONTGENERATION_VALUE_MAX);
    front_trace->gen_table[front_trace->gen_nextfree].trimleft_diff
      = trimleft_diff;
  } else
  {
    front_trace->gen_table[front_trace->gen_nextfree].trimleft_diff = 0;
  }
  front_trace->previoustrimleft = trimleft;
  gt_assert(front_trace->gen_nextfree < front_trace->gen_allocated);
  gt_assert(valid <= GT_FRONTGENERATION_VALUE_MAX);
  front_trace->gen_table[front_trace->gen_nextfree++].valid = valid;
}

void front_trace_add_trace(GtFrontTrace *front_trace,
                           uint8_t backreference,
                           uint32_t localmatch_count)
{
  gt_assert (front_trace != NULL);
  if (front_trace->backref_nextfree >= front_trace->backref_allocated)
  {
    front_trace->backref_allocated
      = front_trace->backref_allocated * 1.2 + 128UL;
    front_trace->backref_table
      = gt_realloc(front_trace->backref_table,
                   sizeof *front_trace->backref_table
                   * front_trace->backref_allocated);
    gt_assert(front_trace->backref_table != NULL);
  }
  gt_assert(front_trace->backref_nextfree < front_trace->backref_allocated);
  front_trace->backref_table[front_trace->backref_nextfree].bits
    = backreference;
  gt_assert(localmatch_count <= front_trace->maxlcs);
  front_trace->backref_table[front_trace->backref_nextfree++].lcs
    = localmatch_count;
#ifdef WITHDISTRIBUTION
  distribution_add(front_trace->lcs_dist,lcs);
#endif
}

/* the following function also works for any point in a front */

static GtUword polished_point2offset(GT_UNUSED const GtFrontTrace *front_trace,
                                     const GtFtPolished_point *pp)
{
  GtWord base_diagonal, pp_diagonal;

  gt_assert(pp != NULL);
  pp_diagonal = (GtWord) pp->alignedlen - (GtWord) GT_MULT2(pp->row);
  gt_assert(pp->distance < front_trace->gen_nextfree);
  base_diagonal = (GtWord) pp->trimleft - (GtWord) pp->distance;
  gt_assert(base_diagonal <= pp_diagonal);
  gt_assert(pp_diagonal < base_diagonal +
                          (GtWord) front_trace->gen_table[pp->distance].valid);
  return (GtUword) (pp_diagonal - base_diagonal);
}

static GtUword valid_total_fronts(const GtFrontGeneration *gen_table,
                                  GtUword start,GtUword end)
{
  GtUword idx, valid_total = 0;

  gt_assert(gen_table);
  for (idx = start; idx < end; idx++)
  {
    valid_total += gen_table[idx].valid;
  }
  return valid_total;
}

static void gt_check_diagonal_run(GT_UNUSED const GtUchar *useq,
                                  GT_UNUSED const GtUchar *vseq,
                                  GT_UNUSED GtWord diagonal,
                                  unsigned int firstrow,
                                  unsigned int nextrow)
{
  GtUword idx;

  gt_assert(useq != NULL && vseq != NULL && firstrow <= nextrow);
  for (idx = firstrow; idx < nextrow; idx++)
  {
    gt_assert (useq[idx] == vseq[idx+diagonal]);
  }
}

static void front_trace2eoplist_directed(GtEoplist *eoplist,
                                         const GtFrontTrace *front_trace,
                                         const GtUchar *useq,
                                         GT_UNUSED GtUword ulen,
                                         const GtUchar *vseq,
                                         GT_UNUSED GtUword vlen,
                                         const GtFtPolished_point *pp)
{
  GtUword distance, localoffset, globaloffset, remainingvalidfronts,
          trimleft;
  GtWord diagonal;
  unsigned int row, lcs;
  uint8_t trace, preferred_eop = FT_EOP_MISMATCH;

  gt_assert(front_trace != NULL && front_trace->gen_nextfree > 0 && pp != NULL);
  localoffset = polished_point2offset(front_trace,pp);
  remainingvalidfronts = valid_total_fronts(front_trace->gen_table,
                                            pp->distance,
                                            front_trace->gen_nextfree);
  gt_assert(remainingvalidfronts <= front_trace->backref_nextfree);
  globaloffset = front_trace->backref_nextfree - remainingvalidfronts;
  distance = pp->distance;
  diagonal = (GtWord) pp->alignedlen - (GtWord) GT_MULT2(pp->row);
  trace = front_trace->backref_table[globaloffset + localoffset].bits;
  lcs = front_trace->backref_table[globaloffset + localoffset].lcs;
  row = pp->row;
  trimleft = pp->trimleft;
  gt_assert(distance < front_trace->gen_nextfree);
  while (distance > 0)
  {
    GtUword nextrowadd;
    GtWord base_diagonal;

    if (eoplist != NULL)
    {
      if (lcs > 0)
      {
        gt_eoplist_match_add(eoplist,lcs);
      }
    } else
    {
      gt_check_diagonal_run(useq, vseq, diagonal, row - lcs, row);
    }
    if (trace & preferred_eop)
    {
      if (preferred_eop == FT_EOP_MISMATCH)
      {
        nextrowadd = 1;
      } else
      {
        if (preferred_eop == FT_EOP_INSERTION)
        {
          gt_assert(-(GtWord) ulen < diagonal);
          diagonal--;
          nextrowadd = 0;
        } else
        {
          gt_assert(preferred_eop == FT_EOP_DELETION);
          gt_assert(diagonal < (GtWord) vlen);
          diagonal++;
          nextrowadd = 1;
        }
      }
    } else
    {
      if (trace & FT_EOP_MISMATCH)
      {
        preferred_eop = FT_EOP_MISMATCH;
        nextrowadd = 1;
      } else
      {
        if (trace & FT_EOP_INSERTION)
        {
          gt_assert(-(GtWord) ulen < diagonal);
          diagonal--;
          preferred_eop = FT_EOP_INSERTION;
          nextrowadd = 0;
        } else
        {
          gt_assert(trace & FT_EOP_DELETION);
          gt_assert(diagonal < (GtWord) vlen);
          diagonal++;
          preferred_eop = FT_EOP_DELETION;
          nextrowadd = 1;
        }
      }
    }
    if (eoplist != NULL)
    {
      if (preferred_eop == FT_EOP_DELETION)
      {
        gt_eoplist_deletion_add(eoplist);
      } else
      {
        if (preferred_eop == FT_EOP_INSERTION)
        {
          gt_eoplist_insertion_add(eoplist);
        } else
        {
          gt_eoplist_mismatch_add(eoplist);
        }
      }
    }
    gt_assert(trimleft >=
              (GtUword) front_trace->gen_table[distance].trimleft_diff);
    trimleft -= (GtUword) front_trace->gen_table[distance].trimleft_diff;
    distance--;
    base_diagonal = (GtWord) trimleft - (GtWord) distance;
    gt_assert(base_diagonal <= diagonal);
    gt_assert(diagonal <
              base_diagonal + (GtWord) front_trace->gen_table[distance].valid);
    localoffset = (GtUword) (diagonal - base_diagonal);
    gt_assert((GtUword) front_trace->gen_table[distance].valid
              <= globaloffset);
    globaloffset -= (GtUword) front_trace->gen_table[distance].valid;
    gt_assert(row >= lcs + nextrowadd);
    row -= lcs + nextrowadd;
    trace = front_trace->backref_table[globaloffset + localoffset].bits;
    lcs = front_trace->backref_table[globaloffset + localoffset].lcs;
  }
  /*printf("avg runlength=%.2f\n",(double) pp->distance/totalrunlength);*/
  gt_assert(globaloffset + localoffset == 0 && trace == 0);
  if (eoplist != NULL && lcs > 0)
  {
    gt_eoplist_match_add(eoplist,lcs);
  }
}

typedef struct
{
  GtUword ulen, vlen;
  GtWord match_score, difference_score;
  bool on_polsize_suffix;
} GtBacktraceFrontInfo;

static GtBacktraceFrontStackelem *stack_top_ptr_get(
                          GtBacktraceFrontStack *stack)
{
  if (stack->nextfree >= stack->allocated)
  {
    stack->allocated = stack->allocated * 1.2 + 128UL;
    stack->space = gt_realloc(stack->space,
                              sizeof *stack->space * stack->allocated);
    gt_assert(stack->space != NULL);
  }
  return stack->space + stack->nextfree++;
}

typedef enum
{
  backtracepath_mismatch,
  backtracepath_deletion,
  backtracepath_insertion
} Backtracepathflag;

static void gt_front_trace_single_push(GtFrontTrace *front_trace,
                                       GtUword match_score,
                                       GtWord diagonal,
                                       GtWord scoresum,
                                       unsigned int row,
                                       GtUword distance,
                                       GtUword globaloffset,
                                       GtUword trimleft,
                                       GtUword lcs_sum,
                                       Backtracepathflag eopcode,
                                       GtUword pathlength)
{
  GtUword localoffset;
  GtWord base_diagonal;
  GtFrontGeneration *gen_table = front_trace->gen_table;
  GtBackreftable *backref_table = front_trace->backref_table;
  GtBacktraceFrontStackelem *stack_top_ptr;

  gt_assert(trimleft >= (GtUword) gen_table[distance+1].trimleft_diff);
  trimleft -= (GtUword) gen_table[distance+1].trimleft_diff;
  base_diagonal = (GtWord) trimleft - (GtWord) distance;
  gt_assert(base_diagonal <= diagonal);
  gt_assert(diagonal < base_diagonal + (GtWord) gen_table[distance].valid);
  localoffset = (GtUword) (diagonal - base_diagonal);
  gt_assert((GtUword) gen_table[distance].valid <= globaloffset);
  globaloffset -= (GtUword) gen_table[distance].valid;
  stack_top_ptr = stack_top_ptr_get(&front_trace->backtracestack);
  stack_top_ptr->diagonal = diagonal;
  stack_top_ptr->distance = distance;
  stack_top_ptr->trace = backref_table[globaloffset + localoffset].bits;
  stack_top_ptr->row = row;
  stack_top_ptr->lcs = backref_table[globaloffset + localoffset].lcs;
  stack_top_ptr->trimleft = trimleft;
  stack_top_ptr->globaloffset = globaloffset;
  stack_top_ptr->lcs_sum = lcs_sum + stack_top_ptr->lcs;
  stack_top_ptr->scoresum = scoresum + stack_top_ptr->lcs * match_score;
  stack_top_ptr->pathlength = pathlength + 1;
  stack_top_ptr->eopcode = eopcode;
}

static void gt_front_trace_backtrace_step(GtBacktraceFrontInfo *bti,
                                          GtFrontTrace *front_trace,
                                          GtWord diagonal,
                                          GtWord scoresum,
                                          GtUword distance,
                                          uint8_t trace,
                                          GtUword globaloffset,
                                          GtUword trimleft,
                                          unsigned int row,
                                          unsigned int lcs,
                                          GtUword lcs_sum,
                                          GtUword pathlength)
{
  gt_assert(distance > 0 && trace != 0);
  if ((trace & FT_EOP_INSERTION) && (!bti->on_polsize_suffix ||
                                     scoresum >= bti->difference_score))
  {
    gt_assert(-(GtWord) bti->ulen < diagonal);
    gt_front_trace_single_push(front_trace,
                               bti->match_score,
                               diagonal - 1,
                               scoresum - bti->difference_score,
                               row - lcs,
                               distance - 1,
                               globaloffset,
                               trimleft,
                               lcs_sum,
                               backtracepath_insertion,
                               pathlength);
    if (!bti->on_polsize_suffix)
    {
      return;
    }
  }
  if ((trace & FT_EOP_DELETION) && (!bti->on_polsize_suffix ||
                                    scoresum >= bti->difference_score))
  {
    gt_assert(diagonal < (GtWord) bti->vlen);
    gt_front_trace_single_push(front_trace,
                               bti->match_score,
                               diagonal + 1,
                               scoresum - bti->difference_score,
                               row - lcs - 1,
                               distance - 1,
                               globaloffset,
                               trimleft,
                               lcs_sum,
                               backtracepath_deletion,
                               pathlength);
    if (!bti->on_polsize_suffix)
    {
      return;
    }
  }
  if ((trace & FT_EOP_MISMATCH) && (!bti->on_polsize_suffix ||
                                    scoresum >= bti->difference_score))
  {
    gt_front_trace_single_push(front_trace,
                               bti->match_score,
                               diagonal,
                               scoresum - bti->difference_score,
                               row - lcs - 1,
                               distance - 1,
                               globaloffset,
                               trimleft,
                               lcs_sum,
                               backtracepath_mismatch,
                               pathlength);
  }
}

static void gt_front_trace_backtracepath2eoplist(GtEoplist *eoplist,
                                                 unsigned int lastlcs,
                                                 const GtBacktraceFrontpath
                                                   *backtracepath,
                                                GtUword elementsinbacktracepath,
                                                GT_UNUSED GtUword ulen,
                                                GT_UNUSED GtUword vlen)
{
  GtUword idx, GT_UNUSED deletions = 0, GT_UNUSED insertions = 0, 
          GT_UNUSED mismatches = 0, GT_UNUSED matches = 0;

  if (lastlcs > 0)
  {
    gt_eoplist_match_add(eoplist,lastlcs);
    matches += lastlcs;
  }
  gt_assert(eoplist != NULL);
  for (idx = 0; idx < elementsinbacktracepath; idx++)
  {
    if (backtracepath[idx].eopcode == backtracepath_deletion)
    {
      gt_eoplist_deletion_add(eoplist);
      deletions++;
    } else
    {
      if (backtracepath[idx].eopcode == backtracepath_insertion)
      {
        gt_eoplist_insertion_add(eoplist);
        insertions++;
      } else
      {
        gt_eoplist_mismatch_add(eoplist);
        mismatches++;
      }
    }
    if (backtracepath[idx].lcs > 0)
    {
      gt_eoplist_match_add(eoplist,backtracepath[idx].lcs);
      matches += backtracepath[idx].lcs;
    }
  }
  /*
  if (matches + mismatches + deletions != ulen)
  {
    fprintf(stderr,
            "matches=" GT_WU ",mismatches=" GT_WU ",deletions=" GT_WU ","
            "sum=" GT_WU " != " GT_WU " = ulen\n",
             matches,mismatches,deletions,
             matches+mismatches+deletions,
             ulen);
  }
  if (matches + mismatches + insertions != vlen)
  {
    fprintf(stderr,
            "matches=" GT_WU ",mismatches=" GT_WU ",insertions=" GT_WU ","
            "sum=" GT_WU " " != " GT_WU " = vlen\n",
             matches,mismatches,insertions,
             matches+mismatches+insertions,
             vlen);
  }
  */
}

static void front_trace2polished_eoplist(GtEoplist *eoplist,
                                         GtFrontTrace *front_trace,
                                         const GtFtPolished_point *pp,
                                         GtUword pol_size,
                                         GtWord match_score,
                                         GtWord difference_score,
                                         const GtUchar *useq,
                                         GtUword ulen,
                                         const GtUchar *vseq,
                                         GtUword vlen)
{
  GtUword localoffset, globaloffset, remainingvalidfronts;
  GtBacktraceFrontStackelem *stack_top_ptr;
  GtBacktraceFrontInfo bti;
  unsigned int lastlcs;

  bti.ulen = ulen;
  bti.vlen = vlen;
  bti.match_score = match_score;
  bti.difference_score = difference_score;
  bti.on_polsize_suffix = true;
  front_trace->backtracestack.nextfree = 0;
  if (front_trace->backtracepath_allocated < pp->distance+1)
  {
    front_trace->backtracepath_allocated = pp->distance + 1;
    front_trace->backtracepath
      = gt_realloc(front_trace->backtracepath,
                   sizeof *front_trace->backtracepath * (pp->distance+1));
  }
  gt_assert(front_trace != NULL && front_trace->gen_nextfree > 0 && pp != NULL);
  localoffset = polished_point2offset(front_trace,pp);
  remainingvalidfronts = valid_total_fronts(front_trace->gen_table,
                                            pp->distance,
                                            front_trace->gen_nextfree);
  gt_assert(remainingvalidfronts <= front_trace->backref_nextfree);
  globaloffset = front_trace->backref_nextfree - remainingvalidfronts;
  stack_top_ptr = stack_top_ptr_get(&front_trace->backtracestack);
  stack_top_ptr->diagonal
    = (GtWord) pp->alignedlen - (GtWord) GT_MULT2(pp->row);
  stack_top_ptr->distance = pp->distance;
  stack_top_ptr->trace
    = front_trace->backref_table[globaloffset + localoffset].bits;
  stack_top_ptr->row = pp->row;
  stack_top_ptr->eopcode = 0;
  lastlcs = stack_top_ptr->lcs
          = front_trace->backref_table[globaloffset + localoffset].lcs;
  stack_top_ptr->scoresum = stack_top_ptr->lcs * match_score;
  stack_top_ptr->globaloffset = globaloffset;
  stack_top_ptr->trimleft = pp->trimleft;
  stack_top_ptr->lcs_sum = stack_top_ptr->lcs;
  stack_top_ptr->pathlength = 0; /* number of errors */
  while (front_trace->backtracestack.nextfree > 0)
  {
    front_trace->backtracestack.nextfree--;
    stack_top_ptr = front_trace->backtracestack.space +
                    front_trace->backtracestack.nextfree;
    if (bti.on_polsize_suffix &&
        stack_top_ptr->lcs_sum + stack_top_ptr->pathlength >= pol_size)
    {
      bti.on_polsize_suffix = false;
    }
    if (stack_top_ptr->pathlength > 0)
    {
      gt_assert(stack_top_ptr->pathlength - 1 <= pp->distance);
      front_trace->backtracepath[stack_top_ptr->pathlength-1].eopcode
        = stack_top_ptr->eopcode;
      front_trace->backtracepath[stack_top_ptr->pathlength-1].lcs
        = stack_top_ptr->lcs;
    }
    if (stack_top_ptr->trace != 0)
    {
      if (eoplist == NULL)
      {
        gt_check_diagonal_run(useq,
                              vseq,
                              stack_top_ptr->diagonal,
                              stack_top_ptr->row -
                              stack_top_ptr->lcs,
                              stack_top_ptr->row);
      }
      gt_front_trace_backtrace_step(&bti,
                                    front_trace,
                                    stack_top_ptr->diagonal,
                                    stack_top_ptr->scoresum,
                                    stack_top_ptr->distance,
                                    stack_top_ptr->trace,
                                    stack_top_ptr->globaloffset,
                                    stack_top_ptr->trimleft,
                                    stack_top_ptr->row,
                                    stack_top_ptr->lcs,
                                    stack_top_ptr->lcs_sum,
                                    stack_top_ptr->pathlength);
    } else
    { /* trace == 0 */
      break;
    }
  }
  gt_assert(stack_top_ptr != NULL);
  if (eoplist != NULL)
  {
    gt_front_trace_backtracepath2eoplist(eoplist,
                                         lastlcs,
                                         front_trace->backtracepath,
                                         stack_top_ptr->pathlength,
                                         ulen,
                                         vlen);
  }
}

void front_trace2eoplist(bool polished,
                         GtEoplist *eoplist,
                         GtFrontTrace *front_trace,
                         const GtFtPolished_point *pp,
                         GtUword pol_size,
                         GtWord match_score,
                         GtWord difference_score,
                         const GtUchar *useq,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vlen)
{
  if (polished)
  {
    front_trace2polished_eoplist(eoplist,
                                 front_trace,
                                 pp,
                                 pol_size,
                                 match_score,
                                 difference_score,
                                 useq,
                                 ulen,
                                 vseq,
                                 vlen);
  } else
  {
    front_trace2eoplist_directed(eoplist,
                                 front_trace,
                                 useq,
                                 ulen,
                                 vseq,
                                 vlen,
                                 pp);
  }
}
void gt_front_trace2eoplist_full_front_directed(GtEoplist *eoplist,
                                                const GtFrontTrace *front_trace,
                                                GtUword distance,
                                                const GtUchar *useq,
                                                GtUword ulen,
                                                const GtUchar *vseq,
                                                GtUword vlen)
{
  const GtBackreftable *basefront, *current;
  GtUword firstindex;
  GtWord diagonal = (GtWord) vlen - (GtWord) ulen;
  uint32_t row;
  uint8_t preferred_eop = FT_EOP_MISMATCH;

  gt_assert(front_trace != NULL &&
            front_trace->backref_nextfree >= 2 * distance + 1);
  basefront = front_trace->backref_table + front_trace->backref_nextfree
                                         - (2 * distance + 1);
  current = basefront + distance + diagonal;
  firstindex = gt_eoplist_length(eoplist);
  gt_assert(ulen <= (GtUword) UINT32_MAX);
  row = ulen;
  while (distance > 0)
  {
    GtUword nextrowadd;

    if (eoplist != NULL)
    {
      if (current->lcs > 0)
      {
        gt_eoplist_match_add(eoplist,current->lcs);
      }
    } else
    {
      gt_check_diagonal_run(useq, vseq, diagonal, row - current->lcs, row);
    }
    if (current->bits & preferred_eop)
    {
      if (preferred_eop == FT_EOP_MISMATCH)
      {
        nextrowadd = 1;
      } else
      {
        if (preferred_eop == FT_EOP_INSERTION)
        {
          gt_assert(-(GtWord) ulen < diagonal);
          diagonal--;
          nextrowadd = 0;
        } else
        {
          gt_assert(preferred_eop == FT_EOP_DELETION);
          gt_assert(diagonal < (GtWord) vlen);
          diagonal++;
          nextrowadd = 1;
        }
      }
    } else
    {
      if (current->bits & FT_EOP_MISMATCH)
      {
        preferred_eop = FT_EOP_MISMATCH;
        nextrowadd = 1;
      } else
      {
        if (current->bits & FT_EOP_INSERTION)
        {
          gt_assert(-(GtWord) ulen < diagonal);
          diagonal--;
          preferred_eop = FT_EOP_INSERTION;
          nextrowadd = 0;
        } else
        {
          gt_assert(current->bits & FT_EOP_DELETION);
          gt_assert(diagonal < (GtWord) vlen);
          diagonal++;
          preferred_eop = FT_EOP_DELETION;
          nextrowadd = 1;
        }
      }
    }
    if (eoplist != NULL)
    {
      if (preferred_eop == FT_EOP_DELETION)
      {
        gt_eoplist_deletion_add(eoplist);
      } else
      {
        if (preferred_eop == FT_EOP_INSERTION)
        {
          gt_eoplist_insertion_add(eoplist);
        } else
        {
          gt_eoplist_mismatch_add(eoplist);
        }
      }
    }
    distance--;
    basefront -= (2 * distance + 1);
    gt_assert(basefront >= front_trace->backref_table);
    gt_assert(row >= current->lcs + nextrowadd);
    row -= current->lcs + nextrowadd;
    current = basefront + distance + diagonal;
  }
  gt_assert(basefront == front_trace->backref_table && current->bits == 0);
  if (eoplist != NULL)
  {
    if (current->lcs > 0)
    {
      gt_eoplist_match_add(eoplist,current->lcs);
    }
    gt_eoplist_reverse_end(eoplist,firstindex);
  }
}
