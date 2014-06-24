/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd

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

#ifndef SEQUENCE_NODE_OUT_VISITOR_H
#define SEQUENCE_NODE_OUT_VISITOR_H

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtSequenceNodeOutVisitor GtSequenceNodeOutVisitor;

#include "core/file_api.h"
#include "extended/genome_node_api.h"
#include "extended/node_visitor_api.h"

GtNodeVisitor* gt_sequence_node_out_visitor_new(GtFile *outfile);
void           gt_sequence_node_out_visitor_keep_sequence_nodes(
                                                     GtSequenceNodeOutVisitor*);
GtUword        gt_sequence_node_out_visitor_node_buffer_size(GtNodeVisitor*);
GtGenomeNode*  gt_sequence_node_out_visitor_get_node(GtNodeVisitor*);
#endif
