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

#ifndef TYPE_NODE_H
#define TYPE_NODE_H

#include "core/bool_matrix.h"

typedef struct GtTypeNode GtTypeNode;

GtTypeNode*   gt_type_node_new(unsigned long num, const char *id);
void          gt_type_node_delete(GtTypeNode*);
unsigned long gt_type_node_num(const GtTypeNode*);
void          gt_type_node_is_a_add(GtTypeNode*, const char*);
const char*   gt_type_node_is_a_get(const GtTypeNode*, unsigned long);
unsigned long gt_type_node_is_a_size(const GtTypeNode*);
void          gt_type_node_part_of_add(GtTypeNode*, const char*);
const char*   gt_type_node_part_of_get(const GtTypeNode*, unsigned long);
unsigned long gt_type_node_part_of_size(const GtTypeNode*);
void          gt_type_node_add_is_a_vertex(GtTypeNode *src,
                                           const GtTypeNode *dst);
bool          gt_type_node_has_parent(GtTypeNode*, const char *id,
                                      GtBoolMatrix *part_of_out_edges,
                                      GtBoolMatrix *part_of_in_edges,
                                      GtArray *node_list);

#endif
