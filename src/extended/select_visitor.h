/*
  Copyright (c) 2007-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef SELECT_VISITOR_H
#define SELECT_VISITOR_H

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtSelectVisitor GtSelectVisitor;

#include "extended/node_visitor.h"
#include "extended/select_stream_api.h"

const GtNodeVisitorClass* gt_select_visitor_class(void);
/* If <strand> is != NUM_OF_GT_STRAND_TYPES, then each genome feature must have
   strand <strand>. */
GtNodeVisitor* gt_select_visitor_new(GtStr *seqid,
                                     GtStr *source,
                                     const GtRange *contain_range,
                                     const GtRange *overlap_range,
                                     GtStrand strand,
                                     GtStrand targetstrand,
                                     bool has_CDS,
                                     unsigned long max_gene_length,
                                     unsigned long max_gene_num,
                                     double min_gene_score,
                                     double max_gene_score,
                                     double min_average_splice_site_prob,
                                     unsigned long feature_num,
                                     GtStrArray *select_files,
                                     GtStr *select_logic,
                                     GtError *err);
void           gt_select_visitor_set_single_intron_factor(GtNodeVisitor*,
                                                          double);
unsigned long  gt_select_visitor_node_buffer_size(GtNodeVisitor*);
GtGenomeNode*  gt_select_visitor_get_node(GtNodeVisitor*);
void           gt_select_visitor_set_drophandler(GtSelectVisitor *fv,
                                                 GtSelectNodeFunc fp,
                                                 void *data);
#endif
