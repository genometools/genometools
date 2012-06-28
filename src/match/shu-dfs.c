/*
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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
#include "core/format64.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/safearith.h"
#include "core/stack-inlined.h"
#include "core/unused_api.h"

#include "match/eis-voiditf.h"
#include "match/shu_unitfile.h"

#include "match/shu-dfs.h"

static inline void add_filenum_count(unsigned long lower,
                                     unsigned long upper,
                                     ShuNode *parent,
                                     BwtSeqpositionextractor *pos_extractor,
                                     unsigned long total_length,
                                     const GtShuUnitFileInfo *unit_info,
                                     const GtEncseq *encseq)
{
  long unsigned row;

  for (row = lower; row < upper; row++)
  {
    unsigned long filenum, position;

    position = gt_BwtSeqpositionextractor_extract(pos_extractor, row);
    if (total_length <= position)
    {
      filenum = unit_info->num_of_files - 1;
    }
    else
    {
      gt_assert((position + 1) <= total_length);
      position = total_length - (position + 1);
      filenum = gt_encseq_filenum(encseq, position);
    }
    /*mapping to genome_num*/
    if (unit_info->map_files != NULL)
    {
      parent->countTermSubtree[0][unit_info->map_files[filenum]]++;
    }
    else
    {
      parent->countTermSubtree[0][filenum]++;
    }
  }
}

static inline unsigned long **get_special_pos(const FMindex *index,
                                         BwtSeqpositionextractor *pos_extractor,
                                         unsigned long num_of_specials)
{
  unsigned long row_idx, special_idx = 0,
                **special_char_rows_and_pos,
                first_special_row;
  GtUchar cc;

  first_special_row = gt_pck_get_nonspecial_count(index);

  gt_array2dim_calloc(special_char_rows_and_pos,
                      2UL,
                      num_of_specials);

  for (row_idx = 0; row_idx < first_special_row; row_idx ++)
  {
    cc = gt_bwtseqgetsymbol(row_idx, index);
    if (ISSPECIAL(cc))
    {
      special_char_rows_and_pos[0][special_idx] = row_idx;
      special_char_rows_and_pos[1][special_idx] =
        gt_BwtSeqpositionextractor_extract(pos_extractor, row_idx);
      special_idx++;
    }
  }
  gt_assert(special_idx == num_of_specials);
  return special_char_rows_and_pos;
}

static inline unsigned long get_start_idx_binary_search(ShuNode *parent,
                                                  unsigned long **special_pos,
                                                  unsigned long max_idx)
{
  unsigned long start_idx = 0;

  if (parent->lower <= special_pos[0][0])
  {
    start_idx = 0;
  }
  else
  {
    if (special_pos[0][max_idx - 1] < parent->lower)
    {
      start_idx = max_idx;
    }
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
        }
        else
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
  return start_idx;
}

static int visit_shu_children(const FMindex *index,
                              ShuNode *parent,
                              GtStackShuNode *stack,
                              const GtEncseq *encseq,
                              Mbtab *tmpmbtab,
                              BwtSeqpositionextractor *pos_extractor,
                              unsigned long *rangeOccs,
                              unsigned long **special_pos,
                              GT_UNUSED unsigned long numofchars,
                              const GtShuUnitFileInfo *unit_info,
                              unsigned long total_length,
                              unsigned long max_idx,
                              GT_UNUSED GtLogger *logger,
                              GT_UNUSED GtError *err)
{
  unsigned long rangesize, idx, num_of_rows;
  unsigned int offset;

  gt_assert(parent->lower < parent->upper);
  num_of_rows = parent->upper - parent->lower;

  rangesize = gt_bwtrangesplitallwithoutspecial(tmpmbtab,
                                                rangeOccs,
                                                index,
                                                parent->lower,
                                                parent->upper);
  gt_assert(rangesize <= numofchars);

  offset = 1U;
  for (idx = 0; idx < rangesize; idx++)
  {
    gt_assert (tmpmbtab[idx].lowerbound <= tmpmbtab[idx].upperbound);
    gt_assert ((tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound) <=
              num_of_rows);
    num_of_rows -= (tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound);

    if (tmpmbtab[idx].lowerbound != tmpmbtab[idx].upperbound)
    {
      if (tmpmbtab[idx].lowerbound + 1 ==
          tmpmbtab[idx].upperbound)
      { /* we found a leaf on parent */
        add_filenum_count(tmpmbtab[idx].lowerbound,
                          tmpmbtab[idx].upperbound,
                          parent,
                          pos_extractor,
                          total_length,
                          unit_info,
                          encseq);
      }
      else
      {
        if (tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound ==
            parent->upper - parent->lower)
        { /* tmpmbtab[idx] is part of a branch and no actual node */
          parent->lower = tmpmbtab[idx].lowerbound;
          parent->upper = tmpmbtab[idx].upperbound;
          parent->depth++;
          return 0;
        }
        else
        { /* tmpmbtab[idx] is a branch of parent node */
          ShuNode *child = NULL;

          GT_STACK_NEXT_FREE(stack,child);
          if (child->countTermSubtree == NULL)
          {
            gt_array2dim_calloc(child->countTermSubtree,
                                numofchars+1UL,
                                unit_info->num_of_genomes);
          }
          else
          {
            unsigned long y_idx, file_idx;
            for (y_idx = 0; y_idx < numofchars+1UL; y_idx++)
            {
              for (file_idx = 0;
                   file_idx < unit_info->num_of_genomes;
                   file_idx++)
              {
                child->countTermSubtree[y_idx][file_idx] = 0;
              }
            }
          }
          child->process = false;
          child->lower = tmpmbtab[idx].lowerbound;
          child->upper = tmpmbtab[idx].upperbound;
          child->parentOffset = offset;
          child->depth = parent->depth + 1;
          offset++;
        }
      }
    }
  }
  if (parent->depth > 0 && num_of_rows > 0) {
    unsigned long start_idx;
    /* unsigned long max_idx - already determined*/

    gt_assert(num_of_rows <= max_idx + 1);
    gt_assert(parent->lower <= special_pos[0][max_idx]);
    gt_assert(special_pos[0][0] <= parent->upper);

    start_idx = get_start_idx_binary_search(parent,
                                            special_pos,
                                            max_idx);
    for (idx = start_idx;
        special_pos[0][idx] < parent->upper && idx <= max_idx;
        idx++)
    {
      unsigned long filenum,
                    position = special_pos[1][idx];

      gt_assert(1UL <= num_of_rows);

      num_of_rows--;

      if (total_length <= position)
      {
        filenum = unit_info->num_of_files - 1;
      }
      else
      {
        gt_assert((position + 1) <= total_length);
        position = total_length - (position + 1);
        filenum = gt_encseq_filenum(encseq, position);
      }
      /*mapping to genome_num*/
      if (unit_info->map_files != NULL)
      {
        parent->countTermSubtree[0][unit_info->map_files[filenum]]++;
      }
      else
      {
        parent->countTermSubtree[0][filenum]++;
      }
    }
    gt_assert(num_of_rows==0);
  }
  parent->process = true;
  return 0;
}

static int process_shu_node(ShuNode *node,
                            GtStackShuNode *stack,
                            uint64_t **shulen,
                            unsigned long num_of_genomes,
                            unsigned long numofchars,
                            GT_UNUSED GtLogger *logger,
                            GT_UNUSED GtError *err)
{
  int had_err = 0;
  uint64_t old;
  unsigned long idx_i, idx_j, termChild_x_i, idx_char;
  unsigned child_c;
  ShuNode *parent = NULL;

  gt_assert(node->process);

  if (node->parentOffset > 0)
  {
    parent = stack->space + stack->nextfree - node->parentOffset;
  }

  for (idx_i = 0; idx_i < num_of_genomes; idx_i++)
  {
    if (!had_err && node->countTermSubtree[0][idx_i] > 0)
    {
      /* scan term */
      termChild_x_i = node->countTermSubtree[0][idx_i];
      for (idx_char = 1UL; idx_char <= numofchars; idx_char++)
      {
        gt_assert(node->countTermSubtree[idx_char][idx_i] <= termChild_x_i);
        termChild_x_i -= node->countTermSubtree[idx_char][idx_i];
      }
      if (termChild_x_i > 0)
      {
        for (idx_j = 0; idx_j < num_of_genomes; idx_j++)
        {
          if (node->countTermSubtree[0][idx_j] > 0 &&
                 idx_j != idx_i)
          {
            old = shulen[idx_i][idx_j];
            shulen[idx_i][idx_j] += ((node->depth + 1) * termChild_x_i);
            if (shulen[idx_i][idx_j] < old)
            {
              had_err = -1;
              gt_error_set(err, "overflow in addition of shuSums! "
                                Formatuint64_t "+ %lu ="
                                Formatuint64_t "\n",
                  PRINTuint64_tcast(old),
                  (node->depth + 1) * termChild_x_i,
                  PRINTuint64_tcast(shulen[idx_i][idx_j]));
            }
          }
        }
      }
      /* scan branch */
      for (child_c = 1U; (unsigned long) child_c <= numofchars; child_c++)
      {
        if (!had_err && node->countTermSubtree[child_c][idx_i] > 0)
        {
          for (idx_j = 0; idx_j < num_of_genomes; idx_j++)
          {
            /* idx_j elem seqIds[x] \ seqIds[child_c] */
            if (node->countTermSubtree[0][idx_j] > 0 &&
                node->countTermSubtree[child_c][idx_j] == 0)
            {
              old = shulen[idx_i][idx_j];
              shulen[idx_i][idx_j] += ((node->depth + 1) *
                                       node->countTermSubtree[child_c][idx_i]);
              if (shulen[idx_i][idx_j] < old)
              {
                had_err = -1;
                gt_error_set(err, "overflow in addition of shuSums! "
                                  Formatuint64_t "+ %lu ="
                                  Formatuint64_t "\n",
                    PRINTuint64_tcast(old),
                    (node->depth + 1) * termChild_x_i,
                    PRINTuint64_tcast(shulen[idx_i][idx_j]));
              }
            }
          }
        }
      }
      if (node->parentOffset > 0)
      {
        gt_assert(parent && parent->countTermSubtree);
        parent->countTermSubtree[0][idx_i] += node->countTermSubtree[0][idx_i];
        parent->countTermSubtree[node->parentOffset][idx_i] =
                                              node->countTermSubtree[0][idx_i];
      }
    }
  }
  return had_err;
}

static int initialise_node(void *node)
{
  int had_err = 0;
  ShuNode *tobeinitialised;
  tobeinitialised = (ShuNode *) node;
  tobeinitialised->countTermSubtree = NULL;
  return had_err;
}

int gt_pck_calculate_shulen(const FMindex *index,
                            const GtShuUnitFileInfo *unit_info,
                            uint64_t **shulen,
                            unsigned long numofchars,
                            unsigned long total_length,
                            GtTimer *timer,
                            GtLogger *logger,
                            GtError *err)
{
  int had_err = 0;
  GtStackShuNode stack;
  ShuNode *root;
  Mbtab *tmpmbtab;
  const unsigned long resize = 64UL;
  unsigned long *rangeOccs,
                **special_char_rows_and_pos,
                depth_idx,
                processed_nodes,
                max_idx = gt_pck_special_occ_in_nonspecial_intervals(index) - 1;
  BwtSeqpositionextractor *pos_extractor;
  const GtEncseq *encseq = unit_info->encseq;

  gt_assert(max_idx < total_length);
  rangeOccs = gt_calloc((size_t) GT_MULT2(numofchars), sizeof (*rangeOccs));
  tmpmbtab = gt_calloc((size_t) (numofchars + 3), sizeof (*tmpmbtab ));
  GT_STACK_INIT_WITH_INITFUNC(&stack, resize, initialise_node);
  pos_extractor = gt_newBwtSeqpositionextractor(index, total_length + 1);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "obtain special pos", stdout);
  }
  special_char_rows_and_pos = get_special_pos(index,
                                              pos_extractor,
                                              max_idx + 1);
  GT_STACK_NEXT_FREE(&stack,root);
  gt_array2dim_calloc(root->countTermSubtree,
                      numofchars+1UL,
                      unit_info->num_of_genomes);
  root->process = false;
  root->parentOffset = 0;
  root->depth = 0;
  root->lower = 0;
  root->upper = total_length + 1;

  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "traverse virtual tree", stdout);
  }
  processed_nodes = 0;
  while (!had_err && !GT_STACK_ISEMPTY(&stack))
  {
    ShuNode *current;

    gt_assert(stack.nextfree > 0);
    current = stack.space + stack.nextfree -1;
    if (current->process)
    {
      GT_STACK_DECREMENTTOP(&stack);
        had_err = process_shu_node(current,
                                   &stack,
                                   shulen,
                                   unit_info->num_of_genomes,
                                   numofchars,
                                   logger,
                                   err);
      processed_nodes++;
    }
    else
    {
      had_err = visit_shu_children(index,
                                   current,
                                   &stack,
                                   encseq,
                                   tmpmbtab,
                                   pos_extractor,
                                   rangeOccs,
                                   special_char_rows_and_pos,
                                   numofchars,
                                   unit_info,
                                   total_length,
                                   max_idx,
                                   logger,
                                   err);
    }
  }
  gt_logger_log(logger, "max stack depth = %lu", GT_STACK_MAXSIZE(&stack));
  gt_log_log("processed nodes= %lu", processed_nodes);
  for (depth_idx = 0; depth_idx < GT_STACK_MAXSIZE(&stack); depth_idx++)
  {
    gt_array2dim_delete(stack.space[depth_idx].countTermSubtree);
  }
  GT_STACK_DELETE(&stack);
  gt_free(rangeOccs);
  gt_free(tmpmbtab);
  gt_freeBwtSeqpositionextractor(pos_extractor);
  gt_array2dim_delete(special_char_rows_and_pos);
  return had_err;
}
