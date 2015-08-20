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
} FrontGeneration;

struct Fronttrace
{
  FrontGeneration *gen_table;
  struct
  {
    uint32_t bits:BACKTRACEBITS,     /* combination of FT_EOP_REPLACEMENT
                                                       FT_EOP_INSERTION
                                                       FT_EOP_DELETION */
             lcs:(32-BACKTRACEBITS); /* longest common suffix */
  } *backref_table;
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

Fronttrace *front_trace_new(void)
{
  Fronttrace *front_trace = gt_malloc(sizeof *front_trace);

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
static size_t front_trace_size(const Fronttrace *front_trace)
{
  gt_assert(front_trace != NULL);
  return sizeof *front_trace->gen_table * front_trace->gen_nextfree +
         sizeof *front_trace->backref_table * front_trace->backref_nextfree;
}
#endif

void front_trace_reset(Fronttrace *front_trace,
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

void front_trace_delete(Fronttrace *front_trace)
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

void front_trace_add_gen(Fronttrace *front_trace,GtUword trimleft,
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

void front_trace_add_trace(Fronttrace *front_trace,uint8_t backreference,
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

static GtUword polished_point2offset(GT_UNUSED const Fronttrace *front_trace,
                                     const Polished_point *pp)
{
  GtWord base_diagonal, pp_diagonal;

  gt_assert(pp != NULL);
  pp_diagonal = (GtWord) pp->alignedlen - (GtWord) GT_MULT2(pp->row);
  gt_assert(pp->distance < front_trace->gen_nextfree);
  base_diagonal = (GtWord) pp->trimleft - (GtWord) pp->distance;
  gt_assert(base_diagonal <= pp_diagonal);
  gt_assert(pp_diagonal <
         base_diagonal + (GtWord) front_trace->gen_table[pp->distance].valid);
  return (GtUword) (pp_diagonal - base_diagonal);
}

static GtUword valid_total_fronts(const FrontGeneration *gen_table,
                                  GtUword start,GtUword end)
{
  GtUword idx, valid_total = 0;

  for (idx = start; idx < end; idx++)
  {
    valid_total += gen_table[idx].valid;
  }
  return valid_total;
}

#ifndef OUTSIDE_OF_GT
void front_trace_multireplacement(GtArrayuint8_t *eoplist,GtUword repnum)
{
  const GtUword addamount = 128;
  while (true)
  {
    if (repnum <= FT_EOPCODE_MAXREPLACEMENT)
    {
      gt_assert(repnum > 0);
      GT_STOREINARRAY(eoplist,uint8_t,addamount,(uint8_t) (repnum - 1));
      break;
    }
    GT_STOREINARRAY(eoplist,uint8_t,addamount,
                    (uint8_t) (FT_EOPCODE_MAXREPLACEMENT - 1));
    repnum -= FT_EOPCODE_MAXREPLACEMENT;
  }
}

void front_trace2eoplist(GtArrayuint8_t *eoplist,
                         const Fronttrace *front_trace,
                         const Polished_point *pp,
                         GT_UNUSED GtUword ulen,
                         GT_UNUSED GtUword vlen)
{
  GtUword distance, localoffset, globaloffset, remainingvalidfronts,
          totalrunlength = 0, trimleft, addamount = 128;
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

    if (lcs > 0)
    {
      front_trace_multireplacement(eoplist,lcs);
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
    if (preferred_eop == FT_EOP_DELETION)
    {
      GT_STOREINARRAY(eoplist,uint8_t,addamount,FT_EOPCODE_DELETION);
    } else
    {
      if (preferred_eop == FT_EOP_INSERTION)
      {
        GT_STOREINARRAY(eoplist,uint8_t,addamount,FT_EOPCODE_INSERTION);
      } else
      {
        GT_STOREINARRAY(eoplist,uint8_t,addamount,0);
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
}
#endif

#define DEBUG
#ifdef DEBUG
static void check_diagonal_run(GT_UNUSED const GtUchar *useq,
                               GT_UNUSED const GtUchar *vseq,
                               GT_UNUSED GtWord diagonal,
                               unsigned int firstrow,
                               unsigned int nextrow)
{
  GtUword idx;

  gt_assert(firstrow <= nextrow);
  for (idx = firstrow; idx < nextrow; idx++)
  {
    gt_assert (useq[idx] == vseq[idx+diagonal]);
  }
}
#endif

void front_trace_verify(const Fronttrace *front_trace,
                        const Polished_point *pp,
                        GT_UNUSED const GtUchar *useq,
                        GT_UNUSED GtUword ulen,
                        GT_UNUSED const GtUchar *vseq,
                        GT_UNUSED GtUword vlen)
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

#ifdef DEBUG
    check_diagonal_run(useq, vseq, diagonal, row - lcs, row);
#endif
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
}

typedef struct
{
  GtWord diagonal;
  GtUword distance,
          globaloffset,
          trimleft;
  unsigned int row, lcs;
  uint8_t trace;
} Backtracestackelem;

typedef struct
{
  Backtracestackelem *space;
  GtUword nextfree, allocated;
} Backtracestack;

static void stack_resize(Backtracestack *stack)
{
  if (stack->nextfree >= stack->allocated)
  {
    stack->allocated = stack->allocated * 1.2 + 128UL;
    stack->space = gt_realloc(stack->space,
                              sizeof *stack->space * stack->allocated);
    gt_assert(stack->space != NULL);
  }
}

static void single_backtrace_step(Backtracestack *stack,
                                  const Fronttrace *front_trace,
                                  GtWord diagonal,
                                  unsigned int row,
                                  GtUword distance,
                                  GtUword globaloffset,
                                  GtUword trimleft)
{
  GtUword localoffset;
  GtWord base_diagonal;

  gt_assert(trimleft >=
            (GtUword) front_trace->gen_table[distance+1].trimleft_diff);
  trimleft -= (GtUword) front_trace->gen_table[distance+1].trimleft_diff;
  base_diagonal = (GtWord) trimleft - (GtWord) distance;
  gt_assert(base_diagonal <= diagonal);
  gt_assert(diagonal <
            base_diagonal + (GtWord) front_trace->gen_table[distance].valid);
  localoffset = (GtUword) (diagonal - base_diagonal);
  gt_assert((GtUword) front_trace->gen_table[distance].valid
            <= globaloffset);
  globaloffset -= (GtUword) front_trace->gen_table[distance].valid;
  stack_resize(stack);
  stack->space[stack->nextfree].diagonal = diagonal;
  stack->space[stack->nextfree].distance = distance;
  stack->space[stack->nextfree].trace
    = front_trace->backref_table[globaloffset + localoffset].bits;
  stack->space[stack->nextfree].row = row;
  stack->space[stack->nextfree].lcs
    = front_trace->backref_table[globaloffset + localoffset].lcs;
  stack->space[stack->nextfree].trimleft = trimleft;
  stack->space[stack->nextfree++].globaloffset = globaloffset;
}

static void backtrace_step(Backtracestack *stack,
                           const Fronttrace *front_trace,
                           GtWord diagonal,
                           GtUword distance,
                           uint8_t trace,
                           GtUword globaloffset,
                           GtUword trimleft,
                           GT_UNUSED unsigned int row,
                           GT_UNUSED unsigned int lcs,
                           GT_UNUSED const GtUchar *useq,
                           GT_UNUSED GtUword ulen,
                           GT_UNUSED const GtUchar *vseq,
                           GT_UNUSED GtUword vlen)
{
  gt_assert(distance > 0 && stack != NULL && trace != 0);
#ifdef DEBUG
  check_diagonal_run(useq, vseq, diagonal, row - lcs, row);
#endif
  if (trace & FT_EOP_INSERTION)
  {
    gt_assert(-(GtWord) ulen < diagonal);
    single_backtrace_step(stack,
                          front_trace,
                          diagonal-1,
                          row - lcs,
                          distance-1,
                          globaloffset,
                          trimleft);
  }
  if (trace & FT_EOP_DELETION)
  {
    gt_assert(diagonal < (GtWord) vlen);
    single_backtrace_step(stack,
                          front_trace,
                          diagonal+1,
                          row - lcs - 1,
                          distance-1,
                          globaloffset,
                          trimleft);
  }
  if (trace & FT_EOP_REPLACEMENT)
  {
    single_backtrace_step(stack,
                          front_trace,
                          diagonal,
                          row - lcs - 1,
                          distance-1,
                          globaloffset,
                          trimleft);
  }
}

void front_trace_verify_all(const Fronttrace *front_trace,
                            const Polished_point *pp,
                            GT_UNUSED const GtUchar *useq,
                            GT_UNUSED GtUword ulen,
                            GT_UNUSED const GtUchar *vseq,
                            GT_UNUSED GtUword vlen)
{
  GtUword localoffset, globaloffset, remainingvalidfronts,
          countofoptimalalignments = 0;
  Backtracestack stack = {NULL,0,0};
  Backtracestackelem *stack_ptr;

  gt_assert(front_trace != NULL && front_trace->gen_nextfree > 0 && pp != NULL);
  localoffset = polished_point2offset(front_trace,pp);
  remainingvalidfronts = valid_total_fronts(front_trace->gen_table,
                                            pp->distance,
                                            front_trace->gen_nextfree);
  gt_assert(remainingvalidfronts <= front_trace->backref_nextfree);
  globaloffset = front_trace->backref_nextfree - remainingvalidfronts;
  stack_resize(&stack);
  stack_ptr = stack.space + stack.nextfree;
  stack_ptr->diagonal = (GtWord) pp->alignedlen - (GtWord) GT_MULT2(pp->row);
  stack_ptr->distance = pp->distance;
  stack_ptr->trace
    = front_trace->backref_table[globaloffset + localoffset].bits;
  stack_ptr->row = pp->row;
  stack_ptr->lcs = front_trace->backref_table[globaloffset + localoffset].lcs;
  stack_ptr->globaloffset = globaloffset;
  stack_ptr->trimleft = pp->trimleft;
  stack.nextfree++;
  while (stack.nextfree > 0)
  {
    stack.nextfree--;
    stack_ptr = stack.space + stack.nextfree;
    if (stack_ptr->trace != 0)
    {
      backtrace_step(&stack,
                     front_trace,
                     stack_ptr->diagonal,
                     stack_ptr->distance,
                     stack_ptr->trace,
                     stack_ptr->globaloffset,
                     stack_ptr->trimleft,
                     stack_ptr->row,
                     stack_ptr->lcs,
                     useq,
                     ulen,
                     vseq,
                     vlen);
    } else
    {
      countofoptimalalignments++;
      if (countofoptimalalignments == 1000)
      {
        break;
      }
    }
  }
  gt_free(stack.space);
}
