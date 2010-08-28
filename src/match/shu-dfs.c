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
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/stack-inlined.h"
#include "core/unused_api.h"

#include "match/eis-voiditf.h"

#include "match/shu-dfs.h"

static inline void add_filenum_count(ShuNode *child,
                                     ShuNode *parent,
                                     BwtSeqpositionextractor *positext,
                                     unsigned long totallength,
                                     unsigned long numoffiles,
                                     const GtEncseq *encseq)
{
  long unsigned row;

  for (row = child->lower; row < child->upper; row++)
  {
    unsigned long filenum, position;

    position = gt_BwtSeqpositionextractor_extract(positext, row);
    if (totallength <= position)
      filenum = numoffiles - 1;
    else
    {
      position = totallength - (position + 1);
      filenum = gt_encseq_filenum(encseq, position);
    }
    parent->countTermSubtree[0][filenum] += 1;
  }
}

static unsigned long **get_special_pos(const FMindex *index,
                                       BwtSeqpositionextractor *positext)
{
  unsigned long row_idx, special_idx = 0,
                **special_char_rows_and_pos,
                num_of_specials,
                first_special_row;

  num_of_specials = gt_pck_special_occ_in_nonspecial_intervals(index);
  first_special_row = gt_pck_get_nonspecial_count(index);

  gt_array2dim_calloc(special_char_rows_and_pos,
                      2UL,
                      num_of_specials);

  for (row_idx = 0; row_idx < first_special_row; row_idx ++)
  {
    GtUchar cc;
    cc = gt_bwtseqgetsymbol(row_idx, index);
    if (ISSPECIAL(cc))
    {
      special_char_rows_and_pos[0][special_idx] = row_idx;
      special_char_rows_and_pos[1][special_idx] =
        gt_BwtSeqpositionextractor_extract(positext, row_idx);
      special_idx += 1;
    }
  }
  gt_assert(num_of_specials == special_idx);
  return special_char_rows_and_pos;
}

static int visit_shu_children(const FMindex *index,
                              ShuNode *parent,
                              GtStackShuNode *stack,
                              const GtEncseq *encseq,
                              Mbtab *tmpmbtab,
                              BwtSeqpositionextractor *positext,
                              unsigned long *rangeOccs,
                              unsigned long **special_pos,
                              unsigned long numofchars,
                              unsigned long numoffiles,
                              unsigned long totallength,
                              unsigned long *stackdepth,
                              bool calculate,
                              GT_UNUSED GtLogger *logger,
                              GT_UNUSED GtError *err)
{
  unsigned long rangesize, idx, num_of_rows;
  unsigned int offset;
  ShuNode child;

  gt_assert(parent->lower < parent->upper);
  num_of_rows = parent->upper - parent->lower;

  if (parent->depth == 0)
  {
    gt_assert(num_of_rows == totallength+1);
    rangesize = gt_bwtrangesplitallwithspecial(tmpmbtab,
                                               rangeOccs,
                                               index,
                                               parent->lower,
                                               parent->upper);
    /* the 3 comes from the seperator, wildcards and terminator */
    gt_assert(rangesize <= numofchars + 3);
  } else
  {
    rangesize = gt_bwtrangesplitallwithoutspecial(tmpmbtab,
                                                  rangeOccs,
                                                  index,
                                                  parent->lower,
                                                  parent->upper);
    gt_assert(rangesize <= numofchars);
  }

  offset = 0;
  for (idx = 0; idx < rangesize; idx++)
  {
    child.process = false;
    child.lower = tmpmbtab[idx].lowerbound;
    child.upper = tmpmbtab[idx].upperbound;
    child.parentOffset = offset + 1;
    child.depth = parent->depth + 1;

    gt_assert (child.lower <= child.upper);
    gt_assert ((child.upper - child.lower) <= num_of_rows);
    num_of_rows = num_of_rows - (child.upper - child.lower);

    if (child.lower == child.upper)
    { /* do nothing, this is no node, this is a missing char */
      continue;
    }
    if (child.lower + 1 == child.upper || numofchars <= idx)
    { /* we found a leave on parent or the interval consists of
      * wildcards or terminators (only in root-interval)*/
      if (calculate)
        add_filenum_count(&child,
                          parent,
                          positext,
                          totallength,
                          numoffiles,
                          encseq);
      continue;
    }
    if (child.upper - child.lower == parent->upper - parent->lower)
    { /* child is part of a branch and no actual node */
      parent->lower = child.lower;
      parent->upper = child.upper;
      parent->depth = child.depth;
      return 0;
    } else
    { /* child is a branch of parent node */
      gt_array2dim_calloc(child.countTermSubtree, 5UL, numoffiles);
      GT_STACK_PUSH(stack, child);
      *stackdepth += 1;
      offset += 1U;
    }
  }
  if (parent->depth != 0 &&
      num_of_rows != 0 && calculate)
  {
    unsigned long start_idx,
                  max_idx = gt_pck_special_occ_in_nonspecial_intervals(index)
                            - 1;

    gt_assert(num_of_rows <= max_idx + 1);
    gt_assert(parent->lower <= special_pos[0][max_idx]);
    gt_assert(special_pos[0][0] <= parent->upper);
    if (parent->lower <= special_pos[0][0])
      start_idx = 0;
    else
    {
      if (special_pos[0][max_idx - 1] < parent->lower)
        start_idx = max_idx;
      else
      {
        unsigned long left, mid, right, len;
        left = 0;
        right = max_idx - 1;
        while (left<=right)
        {
          len = (unsigned long) (right - left);
          mid = left + GT_DIV2(len);
          if (special_pos[0][mid] < parent->lower)
          {
            if (parent->lower <= special_pos[0][mid+1])
            {
              start_idx = mid + 1;
              break;
            }
            left = mid + 1;
          } else
          {
            if (special_pos[0][mid - 1] < parent->lower)
            {
              start_idx = mid;
              break;
            }
            right = mid - 1;
          }
        }
      }
    }
    for (idx = start_idx;
        special_pos[0][idx] < parent->upper && idx <= max_idx;
        idx++)
    {
      unsigned long filenum,
                    position = special_pos[1][idx];

      gt_assert(1UL <= num_of_rows);
      num_of_rows--;

      gt_assert((position + 1) <= totallength);
      if (totallength <= position)
        filenum = numoffiles - 1;
      else
      {
        position = totallength - (position + 1);
        filenum = gt_encseq_filenum(encseq, position);
      }
      parent->countTermSubtree[0][filenum] += 1;
    }
    gt_assert(num_of_rows==0);
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
    parent = stack->space + stack->nextfree - node->parentOffset;
  }

  for (i = 0; i < numoffiles; i++)
  {
    if (node->countTermSubtree[0][i] == 0)
      continue;

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
        shulen[i][j] += ((node->depth + 1) * termChild_x_i);
      }
    }
    /* scan branch */
    for (y = 1U; y < 5U; y++)
    {
      if (node->countTermSubtree[y][i] == 0)
        continue;

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
  return 0;
}

int gt_pck_calculate_shulen(const FMindex *index,
                            const GtEncseq *encseq,
                            double **shulen,
                            unsigned long numofchars,
                            unsigned long totallength,
                            bool calculate,
                            GtProgressTimer *timer,
                            GtLogger *logger,
                            GtError *err)
{
  int had_err = 0;
  GtStackShuNode stack;
  ShuNode root, *current;
  Mbtab *tmpmbtab;
  unsigned long *rangeOccs, **special_char_rows_and_pos;
  unsigned long resize = 64UL; /* XXX make this softcoded */
  unsigned long numoffiles, numofseq, stackdepth, maxdepth;
  BwtSeqpositionextractor *positext;

  rangeOccs = gt_calloc((size_t) GT_MULT2(numofchars), sizeof (*rangeOccs));
  tmpmbtab = gt_calloc((size_t) (numofchars + 3), sizeof (*tmpmbtab ));
  numoffiles = gt_encseq_num_of_files(encseq);
  numofseq = gt_encseq_num_of_sequences(encseq);
  GT_STACK_INIT(&stack, resize);
  positext = gt_newBwtSeqpositionextractor(index, totallength + 1);
  if (timer != NULL)
  {
    gt_progress_timer_start_new_state(timer,
                                      "obtain special pos",
                                      stdout);
  }
  special_char_rows_and_pos = get_special_pos(index, positext);

  gt_array2dim_calloc(root.countTermSubtree, 5UL, numoffiles);
  root.process = false;
  root.parentOffset = 0;
  root.depth = 0;
  root.lower = 0;
  root.upper = totallength + 1;

  if (timer != NULL)
  {
    gt_progress_timer_start_new_state(timer,
                                      "traverse virtual tree",
                                      stdout);
  }
  GT_STACK_PUSH(&stack, root);
  stackdepth = 1UL;
  maxdepth = 0;
  while (!GT_STACK_ISEMPTY(&stack))
  {
    maxdepth = maxdepth < stackdepth ? stackdepth : maxdepth;
    current = &(stack.space[stack.nextfree - 1]);
    if (current->process)
    {
      GT_STACK_DECREMENTTOP(&stack);
      --stackdepth;
      if (calculate)
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
                                   special_char_rows_and_pos,
                                   numofchars,
                                   numoffiles,
                                   totallength,
                                   &stackdepth,
                                   calculate,
                                   logger,
                                   err);
    }
  }
  gt_logger_log(logger, "max stack depth = %lu", maxdepth);
  GT_STACK_DELETE(&stack);
  gt_free(rangeOccs);
  gt_free(tmpmbtab);
  gt_freeBwtSeqpositionextractor(positext);
  gt_array2dim_delete(special_char_rows_and_pos);
  return had_err;
}
