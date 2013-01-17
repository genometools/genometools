/*
  Copyright (c) 2012-2013 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/bool_matrix.h"
#include "core/cstr_api.h"
#include "core/hashmap_api.h"
#include "core/ma_api.h"
#include "core/symbol_api.h"
#include "extended/type_graph.h"
#include "extended/type_node.h"

#define PART_OF           "part_of"
#define MEMBER_OF         "member_of"
#define INTEGRAL_PART_OF  "integral_part_of"

struct GtTypeGraph {
  GtHashmap *name2id, /* maps from name to SO ID */
            *nodemap; /* maps SO ID to actual node */
  GtArray *nodes;
  GtBoolMatrix *part_of_out_edges,
               *part_of_in_edges;
  bool ready;
};

GtTypeGraph* gt_type_graph_new(void)
{
  GtTypeGraph *type_graph = gt_malloc(sizeof (GtTypeGraph));
  type_graph->name2id = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  type_graph->nodemap = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  type_graph->nodes = gt_array_new(sizeof (GtTypeNode*));
  type_graph->part_of_out_edges = gt_bool_matrix_new();
  type_graph->part_of_in_edges = gt_bool_matrix_new();
  type_graph->ready = false;
  return type_graph;
}

void gt_type_graph_delete(GtTypeGraph *type_graph)
{
  unsigned long i;
  if (!type_graph) return;
  gt_bool_matrix_delete(type_graph->part_of_in_edges);
  gt_bool_matrix_delete(type_graph->part_of_out_edges);
  for (i = 0; i < gt_array_size(type_graph->nodes); i++)
    gt_type_node_delete(*(GtTypeNode**) gt_array_get(type_graph->nodes, i));
  gt_array_delete(type_graph->nodes);
  gt_hashmap_delete(type_graph->nodemap);
  gt_hashmap_delete(type_graph->name2id);
  gt_free(type_graph);
}

void gt_type_graph_add_stanza(GtTypeGraph *type_graph,
                              const GtOBOStanza *stanza)
{
  const char *id_value, *name_value;
  unsigned long i, size;
  GtTypeNode *node;
  GtStr *buf;
  gt_assert(type_graph && stanza && !type_graph->ready);
  gt_assert(gt_obo_stanza_size(stanza, "id") == 1);
  gt_assert(gt_obo_stanza_size(stanza, "name") == 1);
  id_value = gt_symbol(gt_obo_stanza_get_value(stanza, "id", 0));
  name_value = gt_symbol(gt_obo_stanza_get_value(stanza, "name", 0));
  gt_assert(id_value);
  gt_assert(name_value);
  gt_assert(!gt_hashmap_get(type_graph->nodemap, id_value));
  node = gt_type_node_new(gt_array_size(type_graph->nodes), id_value);
  gt_hashmap_add(type_graph->name2id, (char*) name_value, (char*) id_value);
  gt_hashmap_add(type_graph->nodemap, (char*) id_value, node);
  gt_array_add(type_graph->nodes, node);
  buf = gt_str_new();
  /* store is_a entries in node, if necessary */
  if ((size = gt_obo_stanza_size(stanza, "is_a"))) {
    for (i = 0; i < size; i++) {
      const char *id = gt_obo_stanza_get_value(stanza, "is_a", i);
      gt_str_reset(buf);
      gt_str_append_cstr_nt(buf, id, strcspn(id, " \n"));
      gt_type_node_is_a_add(node, gt_symbol(gt_str_get(buf)));
    }
  }
  /* store part_of entries in node, if necessary */
  if ((size = gt_obo_stanza_size(stanza, "relationship"))) {
    for (i = 0; i < size; i++) {
      const char *rel = gt_obo_stanza_get_value(stanza, "relationship", i);
      gt_str_reset(buf);
      /* match part_of */
      if (!strncmp(rel, PART_OF, strlen(PART_OF))) {
        const char *part_of = rel + strlen(PART_OF) + 1;
        gt_str_append_cstr_nt(buf, part_of, strcspn(part_of, " \n"));
        gt_type_node_part_of_add(node, gt_symbol(gt_str_get(buf)));
        continue;
      }
      /* match member_of */
      if (!strncmp(rel, MEMBER_OF, strlen(MEMBER_OF))) {
        const char *member_of = rel + strlen(MEMBER_OF) + 1;
        gt_str_append_cstr_nt(buf, member_of, strcspn(member_of, " \n"));
        gt_type_node_part_of_add(node, gt_symbol(gt_str_get(buf)));
        continue;
      }
      /* match integral_part_of */
      if (!strncmp(rel, INTEGRAL_PART_OF, strlen(INTEGRAL_PART_OF))) {
        const char *integral_part_of = rel + strlen(INTEGRAL_PART_OF) + 1;
        gt_str_append_cstr_nt(buf, integral_part_of,
                              strcspn(integral_part_of, " \n"));
        gt_type_node_part_of_add(node, gt_symbol(gt_str_get(buf)));
      }
    }
  }
  gt_str_delete(buf);
}

static void create_vertices(GtTypeGraph *type_graph)
{
  unsigned long i, j;
  GtTypeNode *parent;
  const char *id;
  gt_assert(type_graph && !type_graph->ready);
  /* iterate over nodes */
  for (i = 0; i < gt_array_size(type_graph->nodes); i++) {
    GtTypeNode *node = *(GtTypeNode**) gt_array_get(type_graph->nodes, i);
    /* process is_a parents */
    for (j = 0; j < gt_type_node_is_a_size(node); j++) {
      id = gt_type_node_is_a_get(node, j);
      parent = gt_hashmap_get(type_graph->nodemap, id);
      gt_assert(parent);
      gt_type_node_add_is_a_vertex(node, parent);
    }
    /* process part_of parents */
    for (j = 0; j < gt_type_node_part_of_size(node); j++) {
      id = gt_type_node_part_of_get(node, j);
      parent = gt_hashmap_get(type_graph->nodemap, id);
      gt_assert(parent);
      gt_bool_matrix_set(type_graph->part_of_out_edges, gt_type_node_num(node),
                         gt_type_node_num(parent), true);
      gt_bool_matrix_set(type_graph->part_of_in_edges, gt_type_node_num(parent),
                         gt_type_node_num(node), true);
    }
  }
}

bool gt_type_graph_is_partof(GtTypeGraph *type_graph, const char *parent_type,
                             const char *child_type)
{
  const char *parent_id, *child_id;
  GtTypeNode *child_node;
  gt_assert(type_graph && parent_type && child_type);
  /* make sure graph is built */
  if (!type_graph->ready) {
    create_vertices(type_graph);
    type_graph->ready = true;
  }
  /* get parent ID, if the type is not mappable to an ID, assume it is the ID */
  if (!(parent_id = gt_hashmap_get(type_graph->name2id, parent_type)))
    parent_id = parent_type;
  /* get child ID, if the type is not mappable to an ID, assmue it is the ID */
  if (!(child_id = gt_hashmap_get(type_graph->name2id, child_type)))
    child_id = child_type;
  /* get child node */
  child_node = gt_hashmap_get(type_graph->nodemap, child_id);
  gt_assert(child_node);
  /* check for parent */
  return gt_type_node_has_parent(child_node, parent_id,
                                 type_graph->part_of_out_edges,
                                 type_graph->part_of_in_edges,
                                 type_graph->nodes);
}
