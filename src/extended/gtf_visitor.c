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

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "core/fptr_api.h"
#include "core/unused_api.h"
#include "core/warning.h"
#include "extended/genome_node.h"
#include "extended/genome_visitor_rep.h"
#include "extended/gff3_output.h"
#include "extended/gtf_visitor.h"

struct GTFVisitor {
  const GenomeVisitor parent_instance;
  unsigned long gene_id,
                transcript_id;
  GT_Array *exon_features,
        *CDS_features;
  GenFile *outfp;
};

#define gtf_visitor_cast(GV)\
        genome_visitor_cast(gtf_visitor_class(), GV)

static void gtf_visitor_free(GenomeVisitor *gv)
{
  GTFVisitor *gtf_visitor = gtf_visitor_cast(gv);
  assert(gtf_visitor);
  gt_array_delete(gtf_visitor->exon_features);
  gt_array_delete(gtf_visitor->CDS_features);
}

static int gtf_visitor_comment(GenomeVisitor *gv, GT_Comment *c, GT_UNUSED GT_Error *err)
{
  GTFVisitor *gtf_visitor;
  gt_error_check(err);
  gtf_visitor = gtf_visitor_cast(gv);
  genfile_xprintf(gtf_visitor->outfp, "#%s\n", gt_comment_get_comment(c));
  return 0;
}

static int save_exon_node(GT_GenomeNode *gn, void *data, GT_UNUSED GT_Error *err)
{
  GTFVisitor *gtf_visitor;
  gt_error_check(err);
  assert(gn && data);
  gtf_visitor = (GTFVisitor*) data;
  if (gt_genome_feature_has_type((GT_GenomeFeature*) gn, gft_exon))
    gt_array_add(gtf_visitor->exon_features, gn);
  else if (gt_genome_feature_has_type((GT_GenomeFeature*) gn, gft_CDS))
    gt_array_add(gtf_visitor->CDS_features, gn);
  return 0;
}

static int gtf_show_transcript(GT_GenomeNode *gn, GTFVisitor *gtf_visitor,
                               GT_Error *err)
{
  GT_GenomeFeature *gf;
  unsigned long i;
  int had_err;
  gt_error_check(err);
  assert(gn && gtf_visitor);
  gt_array_reset(gtf_visitor->exon_features);
  gt_array_reset(gtf_visitor->CDS_features);
  had_err = gt_genome_node_traverse_direct_children(gn, gtf_visitor,
                                                 save_exon_node, err);
  if (gt_array_size(gtf_visitor->exon_features)) {
    /* sort exon features */
    qsort(gt_array_get_space(gtf_visitor->exon_features),
          gt_array_size(gtf_visitor->exon_features), sizeof (GT_GenomeNode*),
          (GT_Compare) gt_genome_node_compare);
    /* show exon features */
    gtf_visitor->transcript_id++;
    for (i = 0; i < gt_array_size(gtf_visitor->exon_features); i++) {
      gf = *(GT_GenomeFeature**) gt_array_get(gtf_visitor->exon_features, i);
      gff3_output_leading(gf, gtf_visitor->outfp);
      genfile_xprintf(gtf_visitor->outfp, "gene_id \"%lu\"; transcript_id "
                      "\"%lu.%lu\";\n", gtf_visitor->gene_id,
                      gtf_visitor->gene_id, gtf_visitor->transcript_id);
    }
  }
  if (gt_array_size(gtf_visitor->CDS_features)) {
    /* sort CDS features */
    qsort(gt_array_get_space(gtf_visitor->CDS_features),
          gt_array_size(gtf_visitor->CDS_features), sizeof (GT_GenomeNode*),
          (GT_Compare) gt_genome_node_compare);
    /* show start_codon feature */
    gf = *(GT_GenomeFeature**) gt_array_get(gtf_visitor->CDS_features, 0);
    /* XXX: to be done */

    /* show CDS features */
    for (i = 0; i < gt_array_size(gtf_visitor->CDS_features); i++) {
      gf = *(GT_GenomeFeature**) gt_array_get(gtf_visitor->CDS_features, i);
      gff3_output_leading(gf, gtf_visitor->outfp);
      genfile_xprintf(gtf_visitor->outfp, "gene_id \"%lu\"; transcript_id "
                      "\"%lu.%lu\";\n", gtf_visitor->gene_id,
                      gtf_visitor->gene_id, gtf_visitor->transcript_id);
    }
    /* XXX: show stop_codon feature and shorten last CDS feature */
  }
  return had_err;
}

static int gtf_show_genome_feature(GT_GenomeNode *gn, void *data, GT_Error *err)
{
  GTFVisitor *gtf_visitor = (GTFVisitor*) data;
  GT_GenomeFeature *gf = (GT_GenomeFeature*) gn;
  int had_err = 0;
  if (gt_genome_feature_has_type(gf, gft_gene)) {
      gtf_visitor->gene_id++;
      gtf_visitor->transcript_id = 0;
      had_err = gtf_show_transcript(gn, gtf_visitor, err);
  }
  else if (gt_genome_feature_has_type(gf, gft_mRNA)) {
    had_err = gtf_show_transcript(gn, gtf_visitor, err);
  }
  else if (!(gt_genome_feature_has_type(gf, gft_CDS) ||
             gt_genome_feature_has_type(gf, gft_exon))) {
      warning("skipping GFF3 feature of type \"%s\" (from line %u in file "
              "\"%s\")",
              gt_genome_feature_type_get_cstr(gt_genome_feature_get_type(gf)),
              gt_genome_node_get_line_number(gn), gt_genome_node_get_filename(gn));
  }
  return had_err;
}

static int gtf_visitor_genome_feature(GenomeVisitor *gv, GT_GenomeFeature *gf,
                                      GT_Error *err)
{
  GTFVisitor *gtf_visitor;
  int had_err;
  gt_error_check(err);
  gtf_visitor = gtf_visitor_cast(gv);
  had_err = gt_genome_node_traverse_children((GT_GenomeNode*) gf, gtf_visitor,
                                          gtf_show_genome_feature, false, err);
  return had_err;
}

const GenomeVisitorClass* gtf_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (GTFVisitor),
                                          gtf_visitor_free,
                                          gtf_visitor_comment,
                                          gtf_visitor_genome_feature,
                                          NULL,
                                          NULL };
  return &gvc;
}

GenomeVisitor* gtf_visitor_new(GenFile *outfp)
{
  GenomeVisitor *gv = genome_visitor_create(gtf_visitor_class());
  GTFVisitor *gtf_visitor = gtf_visitor_cast(gv);
  gtf_visitor->gene_id = 0;
  gtf_visitor->exon_features = gt_array_new(sizeof (GT_GenomeNode*));
  gtf_visitor->CDS_features = gt_array_new(sizeof (GT_GenomeNode*));
  gtf_visitor->outfp = outfp;
  return gv;
}
