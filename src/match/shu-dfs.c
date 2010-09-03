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

static inline void add_filenum_count(unsigned long lower,
                                     unsigned long upper,
                                     ShuNode *parent,
                                     BwtSeqpositionextractor *positext,
                                     unsigned long totallength,
                                     unsigned long numoffiles,
                                     const GtEncseq *encseq)
{
  long unsigned row;

  for (row = lower; row < upper; row++)
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
      gt_assert(special_idx <= num_of_specials);
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
                              GT_UNUSED unsigned long numofchars,
                              unsigned long numoffiles,
                              unsigned long totallength,
                              unsigned long *stackdepth,
                              bool calculate,
                              GT_UNUSED GtLogger *logger,
                              GT_UNUSED GtError *err)
{
  unsigned long rangesize, idx, num_of_rows;
  unsigned int offset;
  ShuNode *child;

  gt_assert(parent->lower < parent->upper);
  num_of_rows = parent->upper - parent->lower;

  rangesize = gt_bwtrangesplitallwithoutspecial(tmpmbtab,
                                                rangeOccs,
                                                index,
                                                parent->lower,
                                                parent->upper);
  gt_assert(rangesize <= numofchars);

  offset = 0;
  for (idx = 0; idx < rangesize; idx++)
  {
    gt_assert (tmpmbtab[idx].lowerbound <= tmpmbtab[idx].upperbound);
    gt_assert ((tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound) <=
              num_of_rows);
    num_of_rows = num_of_rows -
                  (tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound);

    if (tmpmbtab[idx].lowerbound != tmpmbtab[idx].upperbound)
    {
      if (tmpmbtab[idx].lowerbound + 1 ==
          tmpmbtab[idx].upperbound)
      { /* we found a leave on parent */
        if (calculate)
          add_filenum_count(tmpmbtab[idx].lowerbound,
                            tmpmbtab[idx].upperbound,
                            parent,
                            positext,
                            totallength,
                            numoffiles,
                            encseq);
      } else
      {
        if (tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound ==
            parent->upper - parent->lower)
        { /* tmpmbtab[idx] is part of a branch and no actual node */
          parent->lower = tmpmbtab[idx].lowerbound;
          parent->upper = tmpmbtab[idx].upperbound;
          parent->depth += 1;
          return 0;
        } else
        { /* tmpmbtab[idx] is a branch of parent node */
          GT_STACK_NEXT_FREE(stack,child);
          if (child->countTermSubtree == NULL)
          {
            gt_array2dim_calloc(child->countTermSubtree, 5UL, numoffiles);
          } else
          {
            unsigned long y_idx, file_idx;
            for (y_idx = 0; y_idx < 5UL; y_idx++)
            {
              for (file_idx = 0; file_idx < numoffiles; file_idx++)
              {
                child->countTermSubtree[y_idx][file_idx] = 0;
              }
            }
          }
          child->process = false;
          child->lower = tmpmbtab[idx].lowerbound;
          child->upper = tmpmbtab[idx].upperbound;
          child->parentOffset = offset + 1;
          child->depth = parent->depth + 1;
          *stackdepth += 1;
          offset += 1U;
        }
      }
    }
  }
  if (parent->depth > 0 &&
      num_of_rows > 0 &&
      calculate)
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
  }
  if (parent->depth != 0 && calculate)
    gt_assert(num_of_rows==0);
  parent->process = true;
  return 0;
}

static int process_shu_node(ShuNode *node,
                            GtStackShuNode *stack,
                            unsigned long **shulen,
                            unsigned long numoffiles,
                            GT_UNUSED GtLogger *logger,
                            GT_UNUSED GtError *err)
{
  unsigned long idx_i, idx_j, termChild_x_i;
  unsigned child_y;
  ShuNode *parent;

  gt_assert(node->process);

  if (node->parentOffset > 0)
  {
    parent = stack->space + stack->nextfree - node->parentOffset;
  }

  for (idx_i = 0; idx_i < numoffiles; idx_i++)
  {
    if (node->countTermSubtree[0][idx_i] > 0)
    {
      /* scan term */
      termChild_x_i = node->countTermSubtree[0][idx_i] -
                      node->countTermSubtree[1][idx_i] -
                      node->countTermSubtree[2][idx_i] -
                      node->countTermSubtree[3][idx_i] -
                      node->countTermSubtree[4][idx_i];
      if (termChild_x_i > 0)
      {
        for (idx_j = 0; idx_j < numoffiles; idx_j++)
        {
          if (node->countTermSubtree[0][idx_j] > 0 &&
                 idx_j != idx_i)
          {
            shulen[idx_i][idx_j] += ((node->depth + 1) * termChild_x_i);
          }
        }
      }
      /* scan branch */
      for (child_y = 1U; child_y < 5U; child_y++)
      {
        if (node->countTermSubtree[child_y][idx_i] > 0)
        {
        for (idx_j = 0; idx_j < numoffiles; idx_j++)
          {
            /* idx_j elem seqIds[x] \ seqIds[child_y] */
            if (node->countTermSubtree[0][idx_j] > 0 &&
                node->countTermSubtree[child_y][idx_j] == 0)
            {
              shulen[idx_i][idx_j] = shulen[idx_i][idx_j] + ((node->depth + 1)
                                     * node->countTermSubtree[child_y][idx_i]);
            }
          }
        }
      }
      if (node->parentOffset > 0)
      {
        parent->countTermSubtree[0][idx_i] += node->countTermSubtree[0][idx_i];
        parent->countTermSubtree[node->parentOffset][idx_i] =
                                              node->countTermSubtree[0][idx_i];
      }
    }
  }
  return 0;
}

int gt_pck_calculate_shulen(const FMindex *index,
                            const GtEncseq *encseq,
                            unsigned long **shulen,
                            unsigned long numofchars,
                            unsigned long totallength,
                            bool calculate,
                            GtProgressTimer *timer,
                            GtLogger *logger,
                            GtError *err)
{
  int had_err = 0;
  GtStackShuNode stack;
  ShuNode *root, *current;
  Mbtab *tmpmbtab;
  unsigned long *rangeOccs, **special_char_rows_and_pos;
  unsigned long resize = 64UL; /* XXX make this softcoded */
  unsigned long numoffiles, numofseq, stackdepth, maxdepth, depth_idx;
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
  GT_STACK_NEXT_FREE(&stack,root);
  gt_array2dim_calloc(root->countTermSubtree, 5UL, numoffiles);
  root->process = false;
  root->parentOffset = 0;
  root->depth = 0;
  root->lower = 0;
  root->upper = totallength + 1;

  if (timer != NULL)
  {
    gt_progress_timer_start_new_state(timer,
                                      "traverse virtual tree",
                                      stdout);
  }
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
  for (depth_idx = 0; depth_idx < maxdepth; depth_idx++)
  {
    gt_array2dim_delete(stack.space[depth_idx].countTermSubtree);
  }
  GT_STACK_DELETE(&stack);
  gt_free(rangeOccs);
  gt_free(tmpmbtab);
  gt_freeBwtSeqpositionextractor(positext);
  gt_array2dim_delete(special_char_rows_and_pos);
  return had_err;
}
