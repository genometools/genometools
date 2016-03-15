/*
  Copyright (c) 2012-2013, 2015 Gordon Gremme <gordon@gremme.org>

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
#include "core/hashmap_api.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "extended/type_node.h"

struct GtTypeNode {
  GtUword num;
  const char *id;
  GtArray *is_a_list,
          *part_of_list;
  GtArray *is_a_out_edges;
  GtHashmap *cache;
  bool transitive_edges_created;
};

GtTypeNode* gt_type_node_new(GtUword num, const char *id)
{
  GtTypeNode *node = gt_calloc(1, sizeof *node);
  node->num = num;
  node->id = id;
  return node;
}

void  gt_type_node_delete(GtTypeNode *type_node)
{
  if (!type_node) return;
  gt_hashmap_delete(type_node->cache);
  gt_array_delete(type_node->is_a_out_edges);
  gt_array_delete(type_node->part_of_list);
  gt_array_delete(type_node->is_a_list);
  gt_free(type_node);
}

GtUword gt_type_node_num(const GtTypeNode *type_node)
{
  gt_assert(type_node);
  return type_node->num;
}

void gt_type_node_is_a_add(GtTypeNode *type_node, const char *id)
{
  gt_assert(type_node && id);
  if (!type_node->is_a_list)
    type_node->is_a_list = gt_array_new(sizeof (const char*));
  gt_array_add(type_node->is_a_list, id);
}

const char* gt_type_node_is_a_get(const GtTypeNode *type_node,
                                  GtUword idx)
{
  gt_assert(type_node && idx < gt_type_node_is_a_size(type_node));
  return *(const char**) gt_array_get(type_node->is_a_list, idx);
}

GtUword gt_type_node_is_a_size(const GtTypeNode *type_node)
{
  gt_assert(type_node);
  return (type_node->is_a_list) ? gt_array_size(type_node->is_a_list) : 0;
}

void gt_type_node_part_of_add(GtTypeNode *type_node, const char *id)
{
  gt_assert(type_node && id);
  if (!type_node->part_of_list)
    type_node->part_of_list = gt_array_new(sizeof (const char*));
  gt_array_add(type_node->part_of_list, id);
}

const char* gt_type_node_part_of_get(const GtTypeNode *type_node,
                                     GtUword idx)
{
  gt_assert(type_node && idx < gt_type_node_part_of_size(type_node));
  return *(const char**) gt_array_get(type_node->part_of_list, idx);
}

GtUword gt_type_node_part_of_size(const GtTypeNode *type_node)
{
  gt_assert(type_node);
  return (type_node->part_of_list) ? gt_array_size(type_node->part_of_list) : 0;
}

void gt_type_node_add_is_a_vertex(GtTypeNode *src, const GtTypeNode *dst)
{
  gt_assert(src && dst);
  if (!src->is_a_out_edges)
    src->is_a_out_edges = gt_array_new(sizeof (GtTypeNode*));
  gt_array_add(src->is_a_out_edges, dst);
}

/* Recursively creates transitive part_of edges for the given <node>.

   That is, if X is_a Y and Z part_of Y, then Z part_of X.

   The method maintains a <node_stack> containing the initial <node> and all
   its 'is_a' children.

   Example from Sequence Ontology (we are looking at mRNA and exon, see
   http://www.sequenceontology.org/browser/current_svn/term/SO:0000234 and
   http://www.sequenceontology.org/browser/current_svn/term/SO:0000147,
   respectively):

   * mRNA is_a mature_transcript
   * mature_transcript is_a transcript
   * transcript_region part_of transcript

   Therefore (transient edges):

   * transcript_region part_of mRNA
   * transcript_region part_of mature_transcript */
static void create_transitive_part_of_edges(GtTypeNode *node,
                                            GtBoolMatrix *part_of_out_edges,
                                            GtBoolMatrix *part_of_in_edges,
                                            GtArray *node_stack)
{
  GtUword i, j;
  /* For every incoming part_of edge of the current node add a corresponding
     part_of edge to all nodes in the stack. Note that the current node is not
     part of the stack yet, therefore this is skipped on the first invocation.
   */
  if (gt_array_size(node_stack)) {
    for (i  = gt_bool_matrix_get_first_column(part_of_in_edges, node->num);
         i != gt_bool_matrix_get_last_column(part_of_in_edges, node->num);
         i  = gt_bool_matrix_get_next_column(part_of_in_edges, node->num, i)) {
      for (j = 0; j < gt_array_size(node_stack); j++) {
        GtTypeNode *child = *(GtTypeNode**) gt_array_get(node_stack, j);
        gt_bool_matrix_set(part_of_out_edges, i, child->num, true);
        gt_bool_matrix_set(part_of_in_edges, child->num, i, true);
      }
    }
  }
  /* add node to stack */
  gt_array_add(node_stack, node);
  /* call method recursively for all is_a edges of the current node */
  for (i = 0; i < gt_array_size(node->is_a_out_edges); i++) {
    GtTypeNode *parent = *(GtTypeNode**) gt_array_get(node->is_a_out_edges, i);
    create_transitive_part_of_edges(parent, part_of_out_edges, part_of_in_edges,
                                    node_stack);
  }
  /* pop node from stack */
  gt_array_pop(node_stack);
}

bool gt_type_node_has_parent(GtTypeNode *node, GtTypeNode *pnode,
                             GtBoolMatrix *part_of_out_edges,
                             GtBoolMatrix *part_of_in_edges,
                             GtArray *node_list, GtHashmap *id2name,
                             unsigned int indentlevel)
{
  GtTypeNode *parent;
  GtUword i;
  bool *result;
  gt_assert(node && pnode);
  gt_log_log("%*scheck if node %s has parent %s", indentlevel * 2, "",
             (const char*) gt_hashmap_get(id2name, node->id),
             (const char*) gt_hashmap_get(id2name, pnode->id));

  /* try cache */
  if (node->cache) {
    if ((result = gt_hashmap_get(node->cache, pnode->id))) {
      gt_log_log("%*sreturn %s (cache hit)", indentlevel * 2, "",
                 *result ?  "true" : "false");
      return *result;
    }
  }
  else
    node->cache = gt_hashmap_new(GT_HASH_DIRECT, NULL, gt_free_func);
  result = gt_malloc(sizeof (bool));

  /* no cache hit found */
  if (node->id == pnode->id) {
    *result = true;
    gt_hashmap_add(node->cache, (char*) pnode->id, result);
    gt_log_log("%*sreturn true", indentlevel * 2, "");
    return true;
  }
  /* create transitive part_of edges for the parent, if they have not been
     created already */
  if (!pnode->transitive_edges_created) {
    GtArray *node_stack;
    node_stack = gt_array_new(sizeof (GtTypeNode*));
    create_transitive_part_of_edges(pnode, part_of_out_edges, part_of_in_edges,
                                    node_stack);
    gt_assert(!gt_array_size(node_stack));
    gt_array_delete(node_stack);
    pnode->transitive_edges_created = true;
  }
  /* traversal of part_of out edges */
  for (i  = gt_bool_matrix_get_first_column(part_of_out_edges, node->num);
       i != gt_bool_matrix_get_last_column(part_of_out_edges, node->num);
       i  = gt_bool_matrix_get_next_column(part_of_out_edges, node->num, i)) {
    parent = *(GtTypeNode**) gt_array_get(node_list, i);
    if (gt_type_node_has_parent(parent, pnode, part_of_out_edges,
                                part_of_in_edges, node_list, id2name,
                                indentlevel + 1)) {
      *result = true;
      gt_hashmap_add(node->cache, (char*) pnode->id, result);
      gt_log_log("%*sreturn true", indentlevel * 2, "");
      return true;
    }
  }
  /* traversal of is_a out edges */
  for (i = 0; i < gt_array_size(node->is_a_out_edges); i++) {
    parent = *(GtTypeNode**) gt_array_get(node->is_a_out_edges, i);
    if (gt_type_node_has_parent(parent, pnode, part_of_out_edges,
                                part_of_in_edges, node_list, id2name,
                                indentlevel + 1)) {
      *result = true;
      gt_hashmap_add(node->cache, (char*) pnode->id, result);
      gt_log_log("%*sreturn true", indentlevel * 2, "");
      return true;
    }
  }
  /* no result found */
  *result = false;
  gt_hashmap_add(node->cache, (char*) pnode->id, result);
  gt_log_log("%*sreturn false", indentlevel * 2, "");
  return false;
}

bool gt_type_node_is_a(GtTypeNode *child_node, const char *parent_id)
{
  GtUword i;
  gt_assert(child_node && parent_id);

  /* XXX: do path compression */
  if (!strcmp(child_node->id, parent_id))
    return true;
  for (i = 0; i < gt_array_size(child_node->is_a_out_edges); i++) {
    GtTypeNode *node = *(GtTypeNode**) gt_array_get(child_node->is_a_out_edges,
                                                    i);
    if (gt_type_node_is_a(node, parent_id))
      return true;
  }
  return false;
}
