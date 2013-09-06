/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_CONTIGS_GRAPH_H
#define RDJ_CONTIGS_GRAPH_H

#include "core/file_api.h"
#include "core/error_api.h"

/* <GtContigsGraph> is a string graph based on the initial contigs
   (reads merged on unambiguous linear paths). */

typedef struct GtContigsGraph GtContigsGraph;

GtContigsGraph* gt_contigs_graph_new(FILE *cjl_i_fp,
                                     FILE *cjl_o_fp,
                                     FILE *junctions_fp,
                                     FILE *rlt_fp,
                                     FILE *depthinfo_fp,
                                     GtError *err);

void            gt_contigs_graph_simplify(GtContigsGraph *cg,
                                          bool restrict_rm_optionals);

void            gt_contigs_graph_extend_contigs(GtContigsGraph *cg,
                                                bool use_only_internal);

void            gt_contigs_graph_show_dot(GtContigsGraph *cg,
                                          GtFile *outfp);

int             gt_contigs_graph_show_dot_subgraph(GtContigsGraph *cg,
                                                   GtFile *outfp,
                                                   GtUword *cnums,
                                                   GtUword nofcnums,
                                                   GtUword maxdepth,
                                                   GtError *err);

void            gt_contigs_graph_enable_dot_show_deleted(GtContigsGraph *cg);

void            gt_contigs_graph_output_paths(GtContigsGraph *cg,
                                              FILE *fp);

void            gt_contigs_graph_delete(GtContigsGraph *cg);

#endif
