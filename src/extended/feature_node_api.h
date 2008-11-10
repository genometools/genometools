/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef FEATURE_NODE_API_H
#define FEATURE_NODE_API_H

#include "core/str_api.h"
#include "core/strand_api.h"
#include "extended/genome_node_api.h"

/* Implements <GtGenomeNode> interface. */
typedef struct GtFeatureNode GtFeatureNode;

/* Create an new <GtFeatureNode*> on sequence with ID <seqid> and type <type>
   which lies from <start> to <end> on strand <strand>.
   The <GtFeatureNode*> stores a new reference to <seqid>, so make sure you do
   not modify the original <seqid> afterwards.
   <start> and <end> always refer to the forward strand, therefore <start> has
   to be smaller or equal than <end>. */
GtGenomeNode* gt_feature_node_new(GtStr *seqid, const char *type,
                                  unsigned long start, unsigned long end,
                                  GtStrand strand);
/* Add <child> node to <parent> node. <parent> takes ownership of <child>.*/
void          gt_feature_node_add_child(GtFeatureNode *parent,
                                        GtFeatureNode *child);

#endif
