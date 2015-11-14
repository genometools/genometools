#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <limits.h>
#ifdef WITHDISTRIBUTION
#include "distribution.h"
#endif
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/unused_api.h"
#include "core/arraydef.h"
#include "core/ma_api.h"
#include "core/types_api.h"
#include "ft-front-generation.h"

#define BACKTRACEBITS 3

typedef struct
{
  uint8_t trimleft_diff,
          valid;
} GtFrontGeneration;

typedef struct
{
  uint32_t bits:BACKTRACEBITS,     /* combination of FT_EOP_REPLACEMENT
                                                     FT_EOP_INSERTION
                                                     FT_EOP_DELETION */
           lcs:(32-BACKTRACEBITS); /* longest common suffix */
} GtBackreftable;

struct GtFronttrace
{
  GtFrontGeneration *gen_table;
  GtBackreftable *backref_table;
#ifdef WITHDISTRIBUTION
  Distribution *lcs_dist, *valid_dist, *trimleft_diff_dist, *space_per_pos_dist;
#endif
  GtUword backref_nextfree,
          backref_allocated,
          gen_nextfree,
          gen_allocated,
          maxlcs,
          previoustrimleft;
};

GtFronttrace *front_trace_new(void)
{
  GtFronttrace *front_trace = gt_malloc(sizeof *front_trace);

  gt_assert(front_trace != NULL);
  /*printf("sizeof *backref_table=" GT_WU " bytes\n",
          sizeof *front_trace->backref_table);*/
  front_trace->backref_table = NULL;
  front_trace->backref_nextfree = 0;
  front_trace->backref_allocated = 0;
  front_trace->gen_table = NULL;
  front_trace->gen_allocated = 0;
  front_trace->gen_nextfree = 0;
  front_trace->previoustrimleft = 0;
  front_trace->maxlcs = (1 << (32-BACKTRACEBITS)) - 1;
#ifdef WITHDISTRIBUTION
  front_trace->lcs_dist = distribution_new();
  front_trace->valid_dist = distribution_new();
  front_trace->trimleft_diff_dist = distribution_new();
  front_trace->space_per_pos_dist = distribution_new();
#endif
  return front_trace;
}

#ifdef WITHDISTRIBUTION
static size_t front_trace_size(const GtFronttrace *front_trace)
{
  gt_assert(front_trace != NULL);
  return sizeof *front_trace->gen_table * front_trace->gen_nextfree +
         sizeof *front_trace->backref_table * front_trace->backref_nextfree;
}
#endif

void front_trace_reset(GtFronttrace *front_trace,
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

void front_trace_delete(GtFronttrace *front_trace)
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
    gt_free(front_trace);
  }
}

void front_trace_add_gen(GtFronttrace *front_trace,GtUword trimleft,
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
    gt_assert(trimleft_diff <= UINT8_MAX);
    front_trace->gen_table[front_trace->gen_nextfree].trimleft_diff
      = trimleft_diff;
  } else
  {
    front_trace->gen_table[front_trace->gen_nextfree].trimleft_diff = 0;
  }
  front_trace->previoustrimleft = trimleft;
  gt_assert(front_trace->gen_nextfree < front_trace->gen_allocated);
  gt_assert(valid <= UINT8_MAX);
  front_trace->gen_table[front_trace->gen_nextfree++].valid = valid;
}

void front_trace_add_trace(GtFronttrace *front_trace,uint8_t backreference,
                           unsigned int lcs)
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
  gt_assert(lcs <= front_trace->maxlcs);
  front_trace->backref_table[front_trace->backref_nextfree++].lcs = lcs;
#ifdef WITHDISTRIBUTION
  distribution_add(front_trace->lcs_dist,lcs);
#endif
}

/* the following function also works for any point in a front */

static GtUword polished_point2offset(GT_UNUSED const GtFronttrace *front_trace,
                                     const Polished_point *pp)
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

char show_eopcode(GT_UNUSED uint8_t eopcode)
{
#ifndef OUTSIDE_OF_GT
  if (eopcode == FT_EOPCODE_DELETION)
  {
    return 'D';
  } else
  {
    if (eopcode == FT_EOPCODE_INSERTION)
    {
      return 'I';
    } else
    {
      return 'R';
    }
  }
#else
  return '\0';
#endif
}

#ifndef OUTSIDE_OF_GT
void eoplist_show(const GtArrayuint8_t *eoplist)
{
  GtUword idx;

  printf("[");
  for (idx = 0; idx < eoplist->nextfreeuint8_t; idx++)
  {
    if (eoplist->spaceuint8_t[idx] == FT_EOPCODE_DELETION)
    {
      printf("D");
    } else
    {
      if (eoplist->spaceuint8_t[idx] == FT_EOPCODE_INSERTION)
      {
        printf("I");
      } else
      {
        GtUword repnum;
        for (repnum = 0; repnum <= eoplist->spaceuint8_t[idx]; repnum++)
        {
          printf("R");
        }
      }
    }
  }
  printf("]\n");
}

#define GT_EOPLIST_PUSH(EOPLIST,EOP)\
        {\
          const GtUword addamount = (EOPLIST)->allocateduint8_t * 0.2 + 128;\
          GT_STOREINARRAY(EOPLIST,uint8_t,addamount,(uint8_t) (EOP));\
        }

void front_trace_multireplacement(GtArrayuint8_t *eoplist,GtUword repnum)
{
  gt_assert(eoplist != NULL && repnum > 0);
  while (true)
  {
    if (repnum <= FT_EOPCODE_MAXREPLACEMENT)
    {
      gt_assert(repnum > 0);
      GT_EOPLIST_PUSH(eoplist,repnum - 1);
      break;
    }
    GT_EOPLIST_PUSH(eoplist,FT_EOPCODE_MAXREPLACEMENT - 1);
    repnum -= FT_EOPCODE_MAXREPLACEMENT;
  }
}
#endif

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

void front_trace2eoplist(GtArrayuint8_t *eoplist,
                         const GtFronttrace *front_trace,
                         const Polished_point *pp,
                         const GtUchar *useq,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vlen)
{
  GtUword distance, localoffset, globaloffset, remainingvalidfronts,
          totalrunlength = 0, trimleft;
  GtWord diagonal;
  unsigned int row, lcs;
  uint8_t trace, preferred_eop = FT_EOP_REPLACEMENT;

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
#ifndef OUTSIDE_OF_GT
      if (lcs > 0)
      {
        front_trace_multireplacement(eoplist,lcs);
      }
#endif
    } else
    {
      gt_check_diagonal_run(useq, vseq, diagonal, row - lcs, row);
    }
    if (trace & preferred_eop)
    {
      totalrunlength++;
      if (preferred_eop == FT_EOP_REPLACEMENT)
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
      if (trace & FT_EOP_REPLACEMENT)
      {
        preferred_eop = FT_EOP_REPLACEMENT;
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
#ifndef OUTSIDE_OF_GT
    if (eoplist != NULL)
    {
      if (preferred_eop == FT_EOP_DELETION)
      {
        GT_EOPLIST_PUSH(eoplist,FT_EOPCODE_DELETION);
      } else
      {
        if (preferred_eop == FT_EOP_INSERTION)
        {
          GT_EOPLIST_PUSH(eoplist,FT_EOPCODE_INSERTION);
        } else
        {
          GT_EOPLIST_PUSH(eoplist,0);
        }
      }
    }
#endif
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
#ifndef OUTSIDE_OF_GT
  if (eoplist != NULL && lcs > 0)
  {
    front_trace_multireplacement(eoplist,lcs);
  }
#endif
}

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
  uint8_t trace, eopcode;
} Backtracestackelem;

typedef struct
{
  Backtracestackelem *space;
  GtUword nextfree, allocated;
} Backtracestack;

typedef struct
{
  uint8_t eopcode;
  unsigned int lcs;
} Backtracepath;

typedef struct
{
  Backtracestack stack;
  const GtFronttrace *front_trace;
  GtUword ulen, vlen;
  GtWord match_score, difference_score;
  bool on_polsize_suffix;
} Backtraceinfo;

static Backtracestackelem *stack_top_ptr_get(Backtracestack *stack)
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

static void single_backtrace_step(Backtraceinfo *bti,
                                  GtWord diagonal,
                                  GtWord scoresum,
                                  unsigned int row,
                                  GtUword distance,
                                  GtUword globaloffset,
                                  GtUword trimleft,
                                  GtUword lcs_sum,
                                  GtUword pathlength,
                                  uint8_t eopcode)
{
  GtUword localoffset;
  GtWord base_diagonal;
  GtFrontGeneration *gen_table = bti->front_trace->gen_table;
  GtBackreftable *backref_table = bti->front_trace->backref_table;
  Backtracestackelem *stack_top_ptr;

  gt_assert(trimleft >= (GtUword) gen_table[distance+1].trimleft_diff);
  trimleft -= (GtUword) gen_table[distance+1].trimleft_diff;
  base_diagonal = (GtWord) trimleft - (GtWord) distance;
  gt_assert(base_diagonal <= diagonal);
  gt_assert(diagonal < base_diagonal + (GtWord) gen_table[distance].valid);
  localoffset = (GtUword) (diagonal - base_diagonal);
  gt_assert((GtUword) gen_table[distance].valid <= globaloffset);
  globaloffset -= (GtUword) gen_table[distance].valid;
  stack_top_ptr = stack_top_ptr_get(&bti->stack);
  stack_top_ptr->diagonal = diagonal;
  stack_top_ptr->distance = distance;
  stack_top_ptr->trace = backref_table[globaloffset + localoffset].bits;
  stack_top_ptr->row = row;
  stack_top_ptr->lcs = backref_table[globaloffset + localoffset].lcs;
  stack_top_ptr->trimleft = trimleft;
  stack_top_ptr->globaloffset = globaloffset;
  stack_top_ptr->lcs_sum = lcs_sum + stack_top_ptr->lcs;
  stack_top_ptr->scoresum = scoresum + stack_top_ptr->lcs * bti->match_score;
  stack_top_ptr->pathlength = pathlength + 1;
  stack_top_ptr->eopcode = eopcode;
}

static void backtrace_step(Backtraceinfo *bti,
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
  uint8_t eopcode = 0;
  gt_assert(distance > 0 && trace != 0);

  if ((trace & FT_EOP_INSERTION) && (!bti->on_polsize_suffix ||
                                     scoresum >= bti->difference_score))
  {
    gt_assert(-(GtWord) bti->ulen < diagonal);
#ifndef OUTSIDE_OF_GT
    eopcode = FT_EOPCODE_INSERTION;
#endif
    single_backtrace_step(bti,
                          diagonal - 1,
                          scoresum - bti->difference_score,
                          row - lcs,
                          distance - 1,
                          globaloffset,
                          trimleft,
                          lcs_sum,
                          pathlength,
                          eopcode);
    if (!bti->on_polsize_suffix)
    {
      return;
    }
  }
  if ((trace & FT_EOP_DELETION) && (!bti->on_polsize_suffix ||
                                    scoresum >= bti->difference_score))
  {
    gt_assert(diagonal < (GtWord) bti->vlen);
#ifndef OUTSIDE_OF_GT
    eopcode = FT_EOPCODE_DELETION;
#endif
    single_backtrace_step(bti,
                          diagonal + 1,
                          scoresum - bti->difference_score,
                          row - lcs - 1,
                          distance - 1,
                          globaloffset,
                          trimleft,
                          lcs_sum,
                          pathlength,
                          eopcode);
    if (!bti->on_polsize_suffix)
    {
      return;
    }
  }
  if ((trace & FT_EOP_REPLACEMENT) && (!bti->on_polsize_suffix ||
                                       scoresum >= bti->difference_score))
  {
#ifndef OUTSIDE_OF_GT
    eopcode = 0;
#endif
    single_backtrace_step(bti,
                          diagonal,
                          scoresum - bti->difference_score,
                          row - lcs - 1,
                          distance - 1,
                          globaloffset,
                          trimleft,
                          lcs_sum,
                          pathlength,
                          eopcode);
  }
}

#ifndef OUTSIDE_OF_GT
static void backtracepath2eoplist(GtArrayuint8_t *eoplist,
                                  unsigned int lastlcs,
                                  const Backtracepath *backtracepath,
                                  GtUword elementsinbacktracepath,
                                  GT_UNUSED GtUword ulen,
                                  GT_UNUSED GtUword vlen)
{
  GtUword idx, deletions = 0, insertions = 0, mismatches = 0, matches = 0;

  if (lastlcs > 0)
  {
    front_trace_multireplacement(eoplist,lastlcs);
    matches += lastlcs;
  }
  for (idx = 0; idx < elementsinbacktracepath; idx++)
  {
    if (backtracepath[idx].eopcode == FT_EOPCODE_DELETION)
    {
      deletions++;
    } else
    {
      if (backtracepath[idx].eopcode == FT_EOPCODE_INSERTION)
      {
        insertions++;
      } else
      {
        mismatches++;
      }
    }
    GT_EOPLIST_PUSH(eoplist,backtracepath[idx].eopcode);
    if (backtracepath[idx].lcs > 0)
    {
      front_trace_multireplacement(eoplist,backtracepath[idx].lcs);
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
#endif

void front_trace2polished_eoplist(GtArrayuint8_t *eoplist,
                                  const GtFronttrace *front_trace,
                                  const Polished_point *pp,
                                  GtUword pol_size,
                                  GtWord match_score,
                                  GtWord difference_score,
                                  const GtUchar *useq,
                                  GtUword ulen,
                                  const GtUchar *vseq,
                                  GtUword vlen)
{
  GtUword localoffset, globaloffset, remainingvalidfronts;
  Backtracestackelem *stack_top_ptr;
  Backtraceinfo bti;
  Backtracepath *backtracepath;
#ifndef OUTSIDE_OF_GT
  unsigned int lastlcs;
#endif
  GT_UNUSED GtUword elementsinbacktracepath = 0;

  bti.stack.space = NULL;
  bti.stack.nextfree = bti.stack.allocated = 0;
  bti.front_trace = front_trace;
  bti.ulen = ulen;
  bti.vlen = vlen;
  bti.match_score = match_score;
  bti.difference_score = difference_score;
  bti.on_polsize_suffix = true;
  backtracepath = gt_malloc(sizeof *backtracepath * (pp->distance+1));
  gt_assert(front_trace != NULL && front_trace->gen_nextfree > 0 && pp != NULL);
  localoffset = polished_point2offset(front_trace,pp);
  remainingvalidfronts = valid_total_fronts(front_trace->gen_table,
                                            pp->distance,
                                            front_trace->gen_nextfree);
  gt_assert(remainingvalidfronts <= front_trace->backref_nextfree);
  globaloffset = front_trace->backref_nextfree - remainingvalidfronts;
  stack_top_ptr = stack_top_ptr_get(&bti.stack);
  stack_top_ptr->diagonal
    = (GtWord) pp->alignedlen - (GtWord) GT_MULT2(pp->row);
  stack_top_ptr->distance = pp->distance;
  stack_top_ptr->trace
    = front_trace->backref_table[globaloffset + localoffset].bits;
  stack_top_ptr->row = pp->row;
#ifndef OUTSIDE_OF_GT
  lastlcs =
#endif
  stack_top_ptr->lcs
    = front_trace->backref_table[globaloffset + localoffset].lcs;
  stack_top_ptr->scoresum = stack_top_ptr->lcs * match_score;
  stack_top_ptr->globaloffset = globaloffset;
  stack_top_ptr->trimleft = pp->trimleft;
  stack_top_ptr->lcs_sum = stack_top_ptr->lcs;
  stack_top_ptr->pathlength = 0; /* number of errors */
  stack_top_ptr->eopcode = 0;
  while (bti.stack.nextfree > 0)
  {
    bti.stack.nextfree--;
    stack_top_ptr = bti.stack.space + bti.stack.nextfree;
    if (bti.on_polsize_suffix &&
        stack_top_ptr->lcs_sum + stack_top_ptr->pathlength >= pol_size)
    {
      bti.on_polsize_suffix = false;
    }
    if (stack_top_ptr->pathlength > 0)
    {
      gt_assert(stack_top_ptr->pathlength - 1 <= pp->distance);
      backtracepath[stack_top_ptr->pathlength-1].eopcode
        = stack_top_ptr->eopcode;
      backtracepath[stack_top_ptr->pathlength-1].lcs = stack_top_ptr->lcs;
      elementsinbacktracepath = stack_top_ptr->pathlength;
    }
    if (stack_top_ptr->trace != 0)
    {
      if (eoplist == NULL)
      {
        gt_check_diagonal_run(useq, vseq,
                              stack_top_ptr->diagonal,
                              stack_top_ptr->row -
                              stack_top_ptr->lcs,
                              stack_top_ptr->row);
      }
      backtrace_step(&bti,
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
#ifndef OUTSIDE_OF_GT
  backtracepath2eoplist(eoplist,
                        lastlcs,
                        backtracepath,
                        elementsinbacktracepath,
                        ulen,
                        vlen);
#endif
  gt_free(backtracepath);
  gt_free(bti.stack.space);
}
