/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include "libgtext/genome_node_info.h"

#define PARENT_MASK             0x0f
#define PARENT_NO_PARENT        0x00 /* the node has not parent */
#define PARENT_ONE_PARENT       0x01 /* the node has exactly one parent */
#define PARENT_MULTIPLE_PARENTS 0x02 /* the node has more then one parent */

#define TREE_MASK               0xf0
#define TREE_UNDETERMINED       0x00 /* the tree status is undetermined */
#define TREE_YES                0x10 /* the a real tree */
#define TREE_NO                 0x20 /* not a tree (node is a DAG instead) */

void genome_node_info_add_parent(GenomeNodeInfo *info)
{
  assert(info);
  switch (*info & PARENT_MASK) {
    case PARENT_NO_PARENT:
      *info |= PARENT_ONE_PARENT;
      break;
    case PARENT_ONE_PARENT:
      *info = (*info & TREE_MASK) | PARENT_MULTIPLE_PARENTS;
      break;
    case PARENT_MULTIPLE_PARENTS:
      /* nothing to do */
      break;
    default: assert(0);
  }
}

bool genome_node_info_multiple_parents(GenomeNodeInfo *info)
{
  assert(info);
  if ((*info & PARENT_MASK) == PARENT_MULTIPLE_PARENTS)
    return true;
  return false;
}

GenomeNodeTreeStatus genome_node_info_get_tree_status(GenomeNodeInfo *info)
{
  GenomeNodeTreeStatus status = GENOME_NODE_TREE_STATUS_UNDETERMINED;
  assert(info);
  switch (*info & TREE_MASK) {
    case TREE_UNDETERMINED:
      status = GENOME_NODE_TREE_STATUS_UNDETERMINED;
      break;
    case TREE_YES:
      status = GENOME_NODE_IS_TREE;
      break;
    case TREE_NO:
      status = GENOME_NODE_IS_NOT_A_TREE;
      break;
    default: assert(0);
  }
  return status;
}

void genome_node_info_set_tree_status(GenomeNodeInfo *info,
                                      GenomeNodeTreeStatus status)
{
  assert(info);
  switch (status) {
    case GENOME_NODE_TREE_STATUS_UNDETERMINED:
      *info = (*info & PARENT_MASK) | TREE_UNDETERMINED;
      break;
    case GENOME_NODE_IS_TREE:
      *info = (*info & PARENT_MASK) | TREE_YES;
      break;
    case GENOME_NODE_IS_NOT_A_TREE:
      *info = (*info & PARENT_MASK) | TREE_NO;
      break;
  }
}
