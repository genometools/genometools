/*
  Copyright (c) 2012-2013 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2012-2013 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/array_api.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/str_array_api.h"
#include "extended/type_node.h"

struct GtTypeNode {
  const char *id;
  GtStrArray *is_a_list,
             *part_of_list;
  GtArray *parents;
};

GtTypeNode* gt_type_node_new(const char *id)
{
  GtTypeNode *node = gt_calloc(1, sizeof *node);
  node->id = id;
  return node;
}

void  gt_type_node_delete(GtTypeNode *type_node)
{
  if (!type_node) return;
  gt_array_delete(type_node->parents);
  gt_str_array_delete(type_node->part_of_list);
  gt_str_array_delete(type_node->is_a_list);
  gt_free(type_node);
}

void gt_type_node_is_a_add(GtTypeNode *type_node, const char *id,
                           unsigned long length)
{
  gt_assert(type_node && id);
  if (!type_node->is_a_list)
    type_node->is_a_list = gt_str_array_new();
  gt_str_array_add_cstr_nt(type_node->is_a_list, id, length);
}

const char* gt_type_node_is_a_get(const GtTypeNode *type_node,
                                  unsigned long idx)
{
  gt_assert(type_node && idx < gt_type_node_is_a_size(type_node));
  return gt_str_array_get(type_node->is_a_list, idx);
}

unsigned long gt_type_node_is_a_size(const GtTypeNode *type_node)
{
  gt_assert(type_node);
  return (type_node->is_a_list) ? gt_str_array_size(type_node->is_a_list) : 0;
}

void gt_type_node_part_of_add(GtTypeNode *type_node, const char *id,
                              unsigned long length)
{
  gt_assert(type_node && id);
  if (!type_node->part_of_list)
    type_node->part_of_list = gt_str_array_new();
  gt_str_array_add_cstr_nt(type_node->part_of_list, id, length);
}

const char* gt_type_node_part_of_get(const GtTypeNode *type_node,
                                  unsigned long idx)
{
  gt_assert(type_node && idx < gt_type_node_part_of_size(type_node));
  return gt_str_array_get(type_node->part_of_list, idx);
}

unsigned long gt_type_node_part_of_size(const GtTypeNode *type_node)
{
  gt_assert(type_node);
  return (type_node->part_of_list) ? gt_str_array_size(type_node->part_of_list)
                                   : 0;
}

void gt_type_node_add_vertex(GtTypeNode *src, const GtTypeNode *dst)
{
  gt_assert(src && dst);
  if (!src->parents)
    src->parents = gt_array_new(sizeof (GtTypeNode*));
  gt_array_add(src->parents, dst);
}

bool gt_type_node_has_parent(GtTypeNode *type_node, const char *id)
{
  unsigned long i;
  gt_assert(type_node && id);
  gt_log_log("check if node %s has parent %s", type_node->id, id);
  if (!strcmp(type_node->id, id)) {
    gt_log_log("return true");
    return true;
  }
  for (i = 0; i < gt_array_size(type_node->parents); i++) {
    GtTypeNode *parent = *(GtTypeNode**) gt_array_get(type_node->parents, i);
    if (gt_type_node_has_parent(parent, id)) {
      gt_log_log("return true");
      return true;
    }
  }
  return false;
}
