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
#include "core/divmodmul.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/stack-inlined.h"
#include "core/unused_api.h"

#include "match/eis-voiditf.h"

#include "match/pck-count-nodes.h"

static void visit_count_children(const FMindex *index,
                           Nodecount *parent,
                           GtStackNodecount *stack,
                           Mbtab *tmpmbtab,
                           unsigned long *rangeOccs,
                           GT_UNUSED unsigned int numofchars)
{
  unsigned long rangesize, idx, num_special;
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
  num_special = parent->upper - parent->lower;
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
      parent->branching--;
      GT_STACK_PUSH(stack, child);
      num_special = 0UL;
      offset++;
    }
    else
    {
      /* we found a leave on parent*/
      if (child.lower + 1 == child.upper)
      {
        parent->leaves++;
        num_special--;
      }
      else
        /* child is a branch of parent node */
      {
        if (child.lower == child.upper)
        {
          /* do nothing, this is no node, this is a missing char */
        }
        else
        {
          GT_STACK_PUSH(stack, child);
          offset++;
          num_special -= (child.upper - child.lower);
        }
      }
    }
  }
  parent->leaves += num_special;
}

static void process_count_node(GtStackNodecount *stack,
                         Nodecount *current)
{
  Nodecount *parent;

  parent = &(stack->space[stack->nextfree - current->parentOffset]);
  parent->leaves += current->leaves;
  parent->branching += current->branching;
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
      {
        /* XXX change to gt_loger_log */
        gt_log_log("on root:\n %lu branching nodes\n %lu leaves\n",
           current->branching, current->leaves);
      }
      else
      {
        process_count_node(&stack, current);
      }
    }
    else
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
