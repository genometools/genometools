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
#include "libgtcore/fptr.h"
#include "libgtcore/unused.h"
#include "libgtcore/warning.h"
#include "libgtext/genome_node.h"
#include "libgtext/genome_visitor_rep.h"
#include "libgtext/gff3_output.h"
#include "libgtext/gtf_visitor.h"

struct GTFVisitor {
  const GenomeVisitor parent_instance;
  unsigned long gene_id,
                transcript_id;
  Array *exon_features,
        *CDS_features;
  GenFile *outfp;
};

#define gtf_visitor_cast(GV)\
        genome_visitor_cast(gtf_visitor_class(), GV)

static void gtf_visitor_free(GenomeVisitor *gv)
{
  GTFVisitor *gtf_visitor = gtf_visitor_cast(gv);
  assert(gtf_visitor);
  array_delete(gtf_visitor->exon_features);
  array_delete(gtf_visitor->CDS_features);
}

static int gtf_visitor_comment(GenomeVisitor *gv, Comment *c, UNUSED Error *err)
{
  GTFVisitor *gtf_visitor;
  error_check(err);
  gtf_visitor = gtf_visitor_cast(gv);
  genfile_xprintf(gtf_visitor->outfp, "#%s\n", comment_get_comment(c));
  return 0;
}

static int save_exon_node(GenomeNode *gn, void *data, UNUSED Error *err)
{
  GTFVisitor *gtf_visitor;
  GenomeFeatureType gft;
  error_check(err);
  assert(gn && data);
  gtf_visitor = (GTFVisitor*) data;
  gft = genome_feature_get_type((GenomeFeature*) gn);
  if (gft == gft_exon)
    array_add(gtf_visitor->exon_features, gn);
  else if (gft == gft_CDS)
    array_add(gtf_visitor->CDS_features, gn);
  return 0;
}

static int gtf_show_transcript(GenomeNode *gn, GTFVisitor *gtf_visitor,
                               Error *e)
{
  GenomeFeature *gf;
  unsigned long i;
  int had_err;
  error_check(e);
  assert(gn && gtf_visitor);
  array_reset(gtf_visitor->exon_features);
  array_reset(gtf_visitor->CDS_features);
  had_err = genome_node_traverse_direct_children(gn, gtf_visitor,
                                                 save_exon_node, e);
  if (array_size(gtf_visitor->exon_features)) {
    /* sort exon features */
    qsort(array_get_space(gtf_visitor->exon_features),
          array_size(gtf_visitor->exon_features), sizeof (GenomeNode*),
          (Compare) genome_node_compare);
    /* show exon features */
    gtf_visitor->transcript_id++;
    for (i = 0; i < array_size(gtf_visitor->exon_features); i++) {
      gf = *(GenomeFeature**) array_get(gtf_visitor->exon_features, i);
      gff3_output_leading(gf, gtf_visitor->outfp);
      genfile_xprintf(gtf_visitor->outfp, "gene_id \"%lu\"; transcript_id "
                      "\"%lu.%lu\";\n", gtf_visitor->gene_id,
                      gtf_visitor->gene_id, gtf_visitor->transcript_id);
    }
  }
  if (array_size(gtf_visitor->CDS_features)) {
    /* sort CDS features */
    qsort(array_get_space(gtf_visitor->CDS_features),
          array_size(gtf_visitor->CDS_features), sizeof (GenomeNode*),
          (Compare) genome_node_compare);
    /* show start_codon feature */
    gf = *(GenomeFeature**) array_get(gtf_visitor->CDS_features, 0);
    /* XXX: to be done */

    /* show CDS features */
    for (i = 0; i < array_size(gtf_visitor->CDS_features); i++) {
      gf = *(GenomeFeature**) array_get(gtf_visitor->CDS_features, i);
      gff3_output_leading(gf, gtf_visitor->outfp);
      genfile_xprintf(gtf_visitor->outfp, "gene_id \"%lu\"; transcript_id "
                      "\"%lu.%lu\";\n", gtf_visitor->gene_id,
                      gtf_visitor->gene_id, gtf_visitor->transcript_id);
    }
    /* XXX: show stop_codon feature and shorten last CDS feature */
  }
  return had_err;
}

static int gtf_show_genome_feature(GenomeNode *gn, void *data, Error *e)
{
  GTFVisitor *gtf_visitor = (GTFVisitor*) data;
  GenomeFeatureType gft;
  int had_err = 0;
  switch ((gft = genome_feature_get_type((GenomeFeature*) gn))) {
    case gft_gene:
      gtf_visitor->gene_id++;
      gtf_visitor->transcript_id = 0;
      had_err = gtf_show_transcript(gn, gtf_visitor, e);
      break;
    case gft_mRNA:
      had_err = gtf_show_transcript(gn, gtf_visitor, e);
      break;
    case gft_CDS:
    case gft_exon:
      /* nothing do do */
      break;
    default:
      warning("skipping GFF3 feature of type \"%s\" (from line %lu in file "
              "\"%s\")", genome_feature_type_get_cstr(gft),
              genome_node_get_line_number(gn), genome_node_get_filename(gn));
  }
  return had_err;
}

static int gtf_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                      Error *e)
{
  GTFVisitor *gtf_visitor;
  int had_err;
  error_check(e);
  gtf_visitor = gtf_visitor_cast(gv);
  had_err = genome_node_traverse_children((GenomeNode*) gf, gtf_visitor,
                                          gtf_show_genome_feature, false, e);
  return had_err;
}

const GenomeVisitorClass* gtf_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (GTFVisitor),
                                          gtf_visitor_free,
                                          gtf_visitor_comment,
                                          gtf_visitor_genome_feature,
                                          NULL };
  return &gvc;
}

GenomeVisitor* gtf_visitor_new(GenFile *outfp)
{
  GenomeVisitor *gv = genome_visitor_create(gtf_visitor_class());
  GTFVisitor *gtf_visitor = gtf_visitor_cast(gv);
  gtf_visitor->gene_id = 0;
  gtf_visitor->exon_features = array_new(sizeof (GenomeNode*));
  gtf_visitor->CDS_features = array_new(sizeof (GenomeNode*));
  gtf_visitor->outfp = outfp;
  return gv;
}
