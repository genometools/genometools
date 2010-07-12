/*
  Copyright (c) 2010 Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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
#include <stdio.h>

#include "core/array2dim_api.h"
#include "core/divmodmul.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/stack-inlined.h"
#include "core/unused_api.h"

#include "match/eis-voiditf.h"

#include "match/shu-dfs.h"

static void visit_count_children(const FMindex *index,
                           Nodecount *parent,
                           GtStackNodecount *stack,
                           Mbtab *tmpmbtab,
                           unsigned long *rangeOccs,
                           unsigned int numofchars)
{
  unsigned long rangesize, idx, num_spezial;
  unsigned int offset;
  Nodecount child;

  gt_assert(parent->lower < parent->upper);
  rangesize = gt_bwtrangesplitallwithoutspecial(tmpmbtab,
                                             rangeOccs,
                                             index,
                                             parent->lower,
                                             parent->upper);
  gt_assert(rangesize <= (unsigned long) numofchars);

  offset = 0U;
  num_spezial = parent->upper - parent->lower;
  for (idx = 0; idx < rangesize; idx++)
  {
    child.lower = tmpmbtab[idx].lowerbound;
    child.upper = tmpmbtab[idx].upperbound;
    child.leaves = 0UL;
    child.branching = 1UL;
    child.parentOffset = offset + 1U;
    child.visited = false;
    child.on_branch = false;

    /* check if child is part of a branch */
    if (child.upper - child.lower == parent->upper - parent->lower)
    {
      /* XXX do sth with non branching nodes */
      child.on_branch = true;
      parent->branching -= 1UL;
      GT_STACK_PUSH(stack, child);
      num_spezial = 0UL;
      offset += 1U;
    } else
    {
      /* we found a leave on parent*/
      if (child.lower + 1 == child.upper)
      {
        parent->leaves += 1;
        num_spezial -= 1;
      } else
        /* child is a branch of parent node */
      {
        if (child.lower == child.upper)
        {
          /* do nothing, this is no node, this is a missing char */
        } else
        {
          GT_STACK_PUSH(stack, child);
          offset += 1U;
          num_spezial -= (child.upper - child.lower);
        }
      }
    }
  }
  parent->leaves += num_spezial;
}

static void process_count_node(GtStackNodecount *stack,
                         Nodecount *current)
{
  Nodecount *parent;

  parent = &(stack->space[stack->nextfree - current->parentOffset]);
  parent->leaves += current->leaves;
  parent->branching += current->branching;
}

static int visit_shu_children(const FMindex *index,
                              ShuNode *parent,
                              GtStackShuNode *stack,
                              const GtEncseq *encseq,
                              Mbtab *tmpmbtab,
                              BwtSeqpositionextractor *positext,
                              unsigned long *rangeOccs,
                              unsigned long numofchars,
                              unsigned long numoffiles,
                              unsigned long totallength,
                              GT_UNUSED GtLogger *logger,
                              GT_UNUSED GtError *err)
{
  unsigned long rangesize, idx, num_of_spezial;
  unsigned int offset;
  ShuNode child;

  gt_assert(parent->lower < parent->upper);
  num_of_spezial = parent->upper - parent->lower;
  rangesize = gt_bwtrangesplitallwithoutspecial(tmpmbtab,
                                             rangeOccs,
                                             index,
                                             parent->lower,
                                             parent->upper);
  gt_assert(rangesize <= numofchars + numoffiles);

  offset = 0U;
  /* check the standard symbols that might lead to child branching node */
  for (idx = 0; idx < rangesize; idx++)
  {
    child.process = false;
    child.lower = tmpmbtab[idx].lowerbound;
    child.upper = tmpmbtab[idx].upperbound;
    child.parentOffset = offset + 1U;
    child.depth = parent->depth + 1;

    num_of_spezial = num_of_spezial - (child.upper - child.lower);

    /* check if child is part of a branch and no actual node */
    if (child.upper - child.lower == parent->upper - parent->lower)
    {
      parent->lower = child.lower;
      parent->upper = child.upper;
      parent->depth = child.depth;
      return 0;
    } else
    {
      /* we found a leave on parent*/
      if (child.lower + 1 == child.upper)
      {
        unsigned long filenum, position;

        position = gt_BwtSeqpositionextractor_extract(positext, child.lower);
        position = totallength - (position + child.depth);
        filenum = gt_encseq_filenum(encseq, position);
        parent->countTermSubtree[0][filenum] += 1;
      } else
        /* child is a branch of parent node */
      {
        if (child.lower == child.upper)
        {
          /* do nothing, this is no node, this is a missing char */
        } else
        {
          gt_array2dim_calloc(child.countTermSubtree, 5UL, numoffiles);
          GT_STACK_PUSH(stack, child);
          offset += 1U;
        }
      }
    }
  }
  /* all remaining "rows" of the parent interval have a special char and are
   * leaves, get the positions and save them to parent */
  for (idx = parent->upper - num_of_spezial; idx < parent->upper; idx++)
  {
    unsigned long filenum, position;

    position = gt_BwtSeqpositionextractor_extract(positext, idx);
    if (child.depth == 1UL)
      position = totallength - position;
    else
      position = totallength - (position + child.depth);
    if (totallength <= position)
      filenum = numoffiles - 1;
    else
    {
      filenum = gt_encseq_filenum(encseq, position);
    }
    parent->countTermSubtree[0][filenum] += 1;
  }
  parent->process = true;
  return 0;
}

static int process_shu_node(ShuNode *node,
                            GtStackShuNode *stack,
                            double **shulen,
                            unsigned long numoffiles,
                            GT_UNUSED GtLogger *logger,
                            GT_UNUSED GtError *err)
{
  unsigned long i, j, termChild_x_i;
  unsigned y;
  ShuNode *parent;

  gt_assert(node->process);

  if (node->parentOffset != 0)
  {
    parent = &(stack->space[stack->nextfree - node->parentOffset]);
  }

  for (i = 0; i < numoffiles; i++)
  {
    termChild_x_i = node->countTermSubtree[0][i] -
                    node->countTermSubtree[1][i] -
                    node->countTermSubtree[2][i] -
                    node->countTermSubtree[3][i] -
                    node->countTermSubtree[4][i];
    /* scan term */
    if (termChild_x_i != 0)
    {
      for (j = 0; j < numoffiles; j++)
      {
        if (node->countTermSubtree[0][j] == 0 ||
               j == i)
          continue;
        shulen[i][j] = shulen[i][j] + ((node->depth + 1) * termChild_x_i);
      }
    }
    /* scan branch */
    for (y = 1U; y < 5U; y++)
    {
      for (j = 0; j < numoffiles; j++)
      {
        /* j elem seqIds[x] \ seqIds[y] */
        if (node->countTermSubtree[0][j] != 0 &&
            node->countTermSubtree[y][j] == 0)
        {
          shulen[i][j] = shulen[i][j] + ((node->depth + 1)
                                          * node->countTermSubtree[y][i]);
        }
      }
    }
    if (node->parentOffset != 0)
    {
      parent->countTermSubtree[0][i] += node->countTermSubtree[0][i];
      parent->countTermSubtree[node->parentOffset][i] =
        node->countTermSubtree[0][i];
    }
  }
  gt_log_log("tiefe: %lu", node->depth);
  if (gt_log_enabled())
  {
    for (i = 0; i < numoffiles; i++)
      {
        printf("debug ");
        for (j = 0; j < numoffiles; j++)
        {
          printf("%f\t", shulen[i][j]);
        }
        printf("\n");
      }
  }
  return 0;
}

int gt_pck_calculate_shulen(const FMindex *index,
                            const GtEncseq *encseq,
                            double **shulen,
                            unsigned long numofchars,
                            unsigned long totallength,
                            GtLogger *logger,
                            GtError *err)
{
  int had_err = 0;
  GtStackShuNode stack;
  ShuNode root, *current;
  Mbtab *tmpmbtab;
  unsigned long *rangeOccs;
  unsigned long resize = 128UL; /* XXX make this softcoded */
  unsigned long numoffiles, numofseq;
  BwtSeqpositionextractor *positext;

  rangeOccs = gt_malloc(sizeof (*rangeOccs) * GT_MULT2(numofchars));
  tmpmbtab = gt_malloc(sizeof (*tmpmbtab ) * numofchars);
  numoffiles = gt_encseq_num_of_files(encseq);
  numofseq = gt_encseq_num_of_sequences(encseq);
  GT_STACK_INIT(&stack, resize);
  positext = gt_newBwtSeqpositionextractor(index, totallength + 1);

  gt_array2dim_calloc(root.countTermSubtree, 5UL, numoffiles);
  root.process = false;
  root.parentOffset = 0;
  root.depth = 0;
  root.lower = 0;
  root.upper = totallength + 1;

  GT_STACK_PUSH(&stack, root);
  while (!GT_STACK_ISEMPTY(&stack))
  {
    current = &(stack.space[stack.nextfree - 1]);
    if (current->process)
    {
      GT_STACK_DECREMENTTOP(&stack);
      had_err = process_shu_node(current,
                                 &stack,
                                 shulen,
                                 numoffiles,
                                 logger,
                                 err);
      gt_array2dim_delete(current->countTermSubtree)
    } else
    {
      had_err = visit_shu_children(index,
                                   current,
                                   &stack,
                                   encseq,
                                   tmpmbtab,
                                   positext,
                                   rangeOccs,
                                   numofchars,
                                   numoffiles,
                                   totallength,
                                   logger,
                                   err);
    }
  }
  GT_STACK_DELETE(&stack);
  gt_free(rangeOccs);
  gt_free(tmpmbtab);
  gt_freeBwtSeqpositionextractor(positext);
  return had_err;
}

void gt_pck_count_nodes_dfs(const FMindex *index,
                        unsigned long totallength,
                        unsigned int numofchars)
{
  GtStackNodecount stack;
  Nodecount root;
  Nodecount *current;
  Mbtab *tmpmbtab;
  unsigned long *rangeOccs;
  unsigned long resize = 128UL; /* XXX make this user definable, or dependable
                                 * on input data */

  GT_STACK_INIT(&stack, resize);
  rangeOccs = gt_malloc(sizeof (*rangeOccs) * GT_MULT2(numofchars));
  tmpmbtab = gt_malloc(sizeof (*tmpmbtab) * numofchars);

  root.lower = 0UL;
  root.upper = totallength + 1;
  root.leaves = 0UL;
  root.branching = 1UL;
  root.parentOffset = 0U;
  root.visited = false;
  root.on_branch = false;

  GT_STACK_PUSH(&stack, root);

  while (!GT_STACK_ISEMPTY(&stack))
  {
    current = &(stack.space[stack.nextfree -1]);
    if (current->visited)
    {
      current = &(GT_STACK_POP(&stack));
      if GT_STACK_ISEMPTY(&stack)
        /* XXX change to gt_loger_log */
        printf("on root:\n %lu branching nodes\n %lu leaves\n",
           current->branching, current->leaves);
      else
      {
        process_count_node(&stack, current);
      }
    } else
    {
      visit_count_children(index, current, &stack,
                     tmpmbtab, rangeOccs, numofchars);
      current->visited = true;
    }
  }
  gt_free(rangeOccs);
  gt_free(tmpmbtab);
  GT_STACK_DELETE(&stack);
}
