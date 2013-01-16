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
#include "core/cstr_api.h"
#include "core/grep_api.h"
#include "core/hashmap_api.h"
#include "core/ma_api.h"
#include "extended/type_graph.h"
#include "extended/type_node.h"

#define PART_OF           "part_of"
#define MEMBER_OF         "member_of"
#define INTEGRAL_PART_OF  "integral_part_of"

struct GtTypeGraph {
  GtHashmap *name2id, /* maps from name to SO ID */
            *nodemap; /* maps SO ID to actual node */
  GtArray *nodes;
  bool ready;
};

GtTypeGraph* gt_type_graph_new(void)
{
  GtTypeGraph *type_graph = gt_malloc(sizeof (GtTypeGraph));
  type_graph->name2id = gt_hashmap_new(GT_HASH_STRING, gt_free_func, NULL);
  type_graph->nodemap = gt_hashmap_new(GT_HASH_STRING, gt_free_func, NULL);
  type_graph->nodes = gt_array_new(sizeof (GtTypeNode*));
  type_graph->ready = false;
  return type_graph;
}

void gt_type_graph_delete(GtTypeGraph *type_graph)
{
  unsigned long i;
  if (!type_graph) return;
  for (i = 0; i < gt_array_size(type_graph->nodes); i++)
    gt_type_node_delete(*(GtTypeNode**) gt_array_get(type_graph->nodes, i));
  gt_array_delete(type_graph->nodes);
  gt_hashmap_delete(type_graph->nodemap);
  gt_hashmap_delete(type_graph->name2id);
  gt_free(type_graph);
}

static bool so_prefix_matches(const char *type)
{
  bool match;
  int rval;
  gt_assert(type);
  rval = gt_grep(&match, "^SO:", type, NULL);
  gt_assert(!rval); /* should not happen */
  return match;
}

void gt_type_graph_add_stanza(GtTypeGraph *type_graph,
                              const GtOBOStanza *stanza)
{
  const char *id_value, *name_value;
  unsigned long i, size;
  char *id_dup;
  GtTypeNode *node;
  gt_assert(type_graph && stanza && !type_graph->ready);
  gt_assert(gt_obo_stanza_size(stanza, "id") == 1);
  gt_assert(gt_obo_stanza_size(stanza, "name") == 1);
  id_value = gt_obo_stanza_get_value(stanza, "id", 0);
  name_value = gt_obo_stanza_get_value(stanza, "name", 0);
  gt_assert(id_value);
  gt_assert(so_prefix_matches(id_value));
  gt_assert(name_value);
  gt_assert(!so_prefix_matches(name_value));
  gt_assert(!gt_hashmap_get(type_graph->nodemap, id_value));
  id_dup = gt_cstr_dup(id_value);
  node = gt_type_node_new(id_dup);
  gt_hashmap_add(type_graph->name2id, gt_cstr_dup(name_value), id_dup);
  gt_hashmap_add(type_graph->nodemap, id_dup, node);
  gt_array_add(type_graph->nodes, node);
  /* store is_a entries in node, if necessary */
  if ((size = gt_obo_stanza_size(stanza, "is_a"))) {
    for (i = 0; i < size; i++) {
      const char *id = gt_obo_stanza_get_value(stanza, "is_a", i);
      gt_type_node_is_a_add(node, id, strcspn(id, " \n"));
    }
  }
  /* store part_of entries in node, if necessary */
  if ((size = gt_obo_stanza_size(stanza, "relationship"))) {
    for (i = 0; i < size; i++) {
      const char *rel = gt_obo_stanza_get_value(stanza, "relationship", i);
      bool match;
      int rval;
      /* match part_of */
      rval = gt_grep(&match, "^"PART_OF, rel, NULL);
      gt_assert(!rval); /* should not happen */
      if (match) {
        const char *part_of = rel + strlen(PART_OF) + 1;
        gt_assert(so_prefix_matches(part_of));
        gt_type_node_part_of_add(node, part_of, strcspn(part_of, " \n"));
        continue;
      }
      /* match member_of */
      rval = gt_grep(&match, "^"MEMBER_OF, rel, NULL);
      gt_assert(!rval); /* should not happen */
      if (match) {
        const char *member_of = rel + strlen(MEMBER_OF) + 1;
        gt_assert(so_prefix_matches(member_of));
        gt_type_node_part_of_add(node, member_of, strcspn(member_of, " \n"));
        continue;
      }
      /* match integral_part_of */
      rval = gt_grep(&match, "^"INTEGRAL_PART_OF, rel, NULL);
      gt_assert(!rval); /* should not happen */
      if (match) {
        const char *integral_part_of = rel + strlen(INTEGRAL_PART_OF) + 1;
        gt_assert(so_prefix_matches(integral_part_of));
        gt_type_node_part_of_add(node, integral_part_of,
                                 strcspn(integral_part_of, " \n"));
      }
    }
  }
}

static void create_vertices(GtTypeGraph *type_graph)
{
  unsigned long i, j;
  GtTypeNode *parent;
  const char *id;
  gt_assert(type_graph && !type_graph->ready);
  gt_assert(gt_hashmap_get(type_graph->nodemap, "SO:0000110"));
  /* iterate over nodes */
  for (i = 0; i < gt_array_size(type_graph->nodes); i++) {
    GtTypeNode *node = *(GtTypeNode**) gt_array_get(type_graph->nodes, i);
    /* process is_a parents */
    for (j = 0; j < gt_type_node_is_a_size(node); j++) {
      id = gt_type_node_is_a_get(node, j);
      parent = gt_hashmap_get(type_graph->nodemap, id);
      gt_assert(parent);
      gt_type_node_add_vertex(node, parent);
    }
    /* process part_of parents */
    for (j = 0; j < gt_type_node_part_of_size(node); j++) {
      id = gt_type_node_part_of_get(node, j);
      parent = gt_hashmap_get(type_graph->nodemap, id);
      gt_assert(parent);
      gt_type_node_add_vertex(node, parent);
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
  /* get parent ID */
  if (!so_prefix_matches(parent_type)) {
    parent_id = gt_hashmap_get(type_graph->name2id, parent_type);
    gt_assert(parent_id);
  }
  else
    parent_id = parent_type;
  /* get child ID */
  if (!so_prefix_matches(child_type)) {
    child_id = gt_hashmap_get(type_graph->name2id, child_type);
    gt_assert(child_id);
  }
  else
    child_id = child_type;
  /* get child node */
  child_node = gt_hashmap_get(type_graph->nodemap, child_id);
  gt_assert(child_node);
  /* check for parent */
  return gt_type_node_has_parent(child_node, parent_id);
}
