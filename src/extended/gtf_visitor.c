/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/fptr_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_output.h"
#include "extended/gtf_visitor.h"
#include "extended/node_visitor_api.h"

struct GtGTFVisitor {
  const GtNodeVisitor parent_instance;
  unsigned long gene_id,
                transcript_id;
  GtArray *exon_features,
          *CDS_features;
  GtFile *outfp;
};

#define gtf_visitor_cast(GV)\
        gt_node_visitor_cast(gt_gtf_visitor_class(), GV)

static void gtf_visitor_free(GtNodeVisitor *nv)
{
  GtGTFVisitor *gtf_visitor = gtf_visitor_cast(nv);
  gt_assert(gtf_visitor);
  gt_array_delete(gtf_visitor->exon_features);
  gt_array_delete(gtf_visitor->CDS_features);
}

static int gtf_visitor_comment_node(GtNodeVisitor *nv, GtCommentNode *c,
                                    GT_UNUSED GtError *err)
{
  GtGTFVisitor *gtf_visitor;
  gt_error_check(err);
  gtf_visitor = gtf_visitor_cast(nv);
  gt_file_xprintf(gtf_visitor->outfp, "#%s\n", gt_comment_node_get_comment(c));
  return 0;
}

static int save_exon_node(GtFeatureNode *fn, void *data, GT_UNUSED GtError *err)
{
  GtGTFVisitor *gtf_visitor;
  gt_error_check(err);
  gt_assert(fn && data);
  gtf_visitor = (GtGTFVisitor*) data;
  if (gt_feature_node_has_type(fn, gt_ft_exon))
    gt_array_add(gtf_visitor->exon_features, fn);
  else if (gt_feature_node_has_type(fn, gt_ft_CDS))
    gt_array_add(gtf_visitor->CDS_features, fn);
  return 0;
}

static int gtf_show_transcript(GtFeatureNode *feature_node,
                               GtGTFVisitor *gtf_visitor, GtError *err)
{
  GtFeatureNode *fn;
  unsigned long i;
  int had_err;
  gt_error_check(err);
  gt_assert(feature_node && gtf_visitor);
  gt_array_reset(gtf_visitor->exon_features);
  gt_array_reset(gtf_visitor->CDS_features);
  had_err = gt_feature_node_traverse_direct_children(feature_node, gtf_visitor,
                                                     save_exon_node, err);
  if (gt_array_size(gtf_visitor->exon_features)) {
    /* sort exon features */
    qsort(gt_array_get_space(gtf_visitor->exon_features),
          gt_array_size(gtf_visitor->exon_features), sizeof (GtGenomeNode*),
          (GtCompare) gt_genome_node_compare);
    /* show exon features */
    gtf_visitor->transcript_id++;
    for (i = 0; i < gt_array_size(gtf_visitor->exon_features); i++) {
      fn = *(GtFeatureNode**) gt_array_get(gtf_visitor->exon_features, i);
      gt_gff3_output_leading(fn, gtf_visitor->outfp);
      gt_file_xprintf(gtf_visitor->outfp, "gene_id \"%lu\"; transcript_id "
                      "\"%lu.%lu\";\n", gtf_visitor->gene_id,
                      gtf_visitor->gene_id, gtf_visitor->transcript_id);
    }
  }
  if (gt_array_size(gtf_visitor->CDS_features)) {
    /* sort CDS features */
    qsort(gt_array_get_space(gtf_visitor->CDS_features),
          gt_array_size(gtf_visitor->CDS_features), sizeof (GtGenomeNode*),
          (GtCompare) gt_genome_node_compare);
    /* show start_codon feature */
    fn = *(GtFeatureNode**) gt_array_get(gtf_visitor->CDS_features, 0);
    /* XXX: to be done */

    /* show CDS features */
    for (i = 0; i < gt_array_size(gtf_visitor->CDS_features); i++) {
      fn = *(GtFeatureNode**) gt_array_get(gtf_visitor->CDS_features, i);
      gt_gff3_output_leading(fn, gtf_visitor->outfp);
      gt_file_xprintf(gtf_visitor->outfp, "gene_id \"%lu\"; transcript_id "
                      "\"%lu.%lu\";\n", gtf_visitor->gene_id,
                      gtf_visitor->gene_id, gtf_visitor->transcript_id);
    }
    /* XXX: show stop_codon feature and shorten last CDS feature */
  }
  return had_err;
}

static int gtf_show_feature_node(GtFeatureNode *fn, void *data, GtError *err)
{
  GtGTFVisitor *gtf_visitor = (GtGTFVisitor*) data;
  int had_err = 0;
  if (gt_feature_node_has_type(fn, gt_ft_gene)) {
      gtf_visitor->gene_id++;
      gtf_visitor->transcript_id = 0;
      had_err = gtf_show_transcript(fn, gtf_visitor, err);
  }
  else if (gt_feature_node_has_type(fn, gt_ft_mRNA)) {
    had_err = gtf_show_transcript(fn, gtf_visitor, err);
  }
  else if (!(gt_feature_node_has_type(fn, gt_ft_CDS) ||
             gt_feature_node_has_type(fn, gt_ft_exon))) {
      gt_warning("skipping GFF3 feature of type \"%s\" (from line %u in file "
                 "\"%s\")",
                 gt_feature_node_get_type(fn),
                 gt_genome_node_get_line_number((GtGenomeNode*) fn),
                 gt_genome_node_get_filename((GtGenomeNode*) fn));
  }
  return had_err;
}

static int gtf_visitor_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                    GtError *err)
{
  GtGTFVisitor *gtf_visitor;
  int had_err;
  gt_error_check(err);
  gtf_visitor = gtf_visitor_cast(nv);
  had_err = gt_feature_node_traverse_children(fn, gtf_visitor,
                                              gtf_show_feature_node, false,
                                              err);
  return had_err;
}

const GtNodeVisitorClass* gt_gtf_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtGTFVisitor),
                                    gtf_visitor_free,
                                    gtf_visitor_comment_node,
                                    gtf_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_gtf_visitor_new(GtFile *outfp)
{
  GtNodeVisitor *nv = gt_node_visitor_create(gt_gtf_visitor_class());
  GtGTFVisitor *gtf_visitor = gtf_visitor_cast(nv);
  gtf_visitor->gene_id = 0;
  gtf_visitor->exon_features = gt_array_new(sizeof (GtGenomeNode*));
  gtf_visitor->CDS_features = gt_array_new(sizeof (GtGenomeNode*));
  gtf_visitor->outfp = outfp;
  return nv;
}
