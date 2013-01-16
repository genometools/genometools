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

#ifndef TYPE_GRAPH_H
#define TYPE_GRAPH_H

#include "extended/obo_stanza.h"

typedef struct GtTypeGraph GtTypeGraph;

GtTypeGraph* gt_type_graph_new(void);
void         gt_type_graph_delete(GtTypeGraph *type_graph);
void         gt_type_graph_add_stanza(GtTypeGraph *type_graph,
                                      const GtOBOStanza *obo_stanza);
bool         gt_type_graph_is_partof(GtTypeGraph *type_graph,
                                     const char *parent_type,
                                     const char *child_type);

#endif
