/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <libgtext/genome_node.h>
#include <libgtext/genome_visitor_rep.h>
#include <libgtext/gtf_visitor.h>

struct GTFVisitor {
  const GenomeVisitor parent_instance;
  unsigned long gene_id,
                transcript_id;
  Array *exon_features;
  GenFile *outfp;
};

#define gtf_visitor_cast(GV)\
        genome_visitor_cast(gtf_visitor_class(), GV)

static void gtf_visitor_free(GenomeVisitor *gv, Env *env)
{
  GTFVisitor *gtf_visitor = gtf_visitor_cast(gv);
  assert(gtf_visitor);
  array_delete(gtf_visitor->exon_features, env);
}

static int gtf_visitor_comment(GenomeVisitor *gv, Comment *c, Env *env)
{
  GTFVisitor *gtf_visitor;
  env_error_check(env);
  gtf_visitor = gtf_visitor_cast(gv);
  genfile_xprintf(gtf_visitor->outfp, "#%s\n", comment_get_comment(c));
  return 0;
}

static int save_exon_node(GenomeNode *gn, void *data, Env *env)
{
  GTFVisitor *gtf_visitor;
  env_error_check(env);
  assert(gn && data);
  gtf_visitor = (GTFVisitor*) data;
  if (genome_feature_get_type((GenomeFeature*) gn) == gft_exon)
    array_add(gtf_visitor->exon_features, gn, env);
  return 0;
}

static int gtf_show_transcript(GenomeNode *gn, GTFVisitor *gtf_visitor,
                               Env *env)
{
  int has_err;
  env_error_check(env);
  assert(gn && gtf_visitor);
  array_set_size(gtf_visitor->exon_features, 0);
  has_err = genome_node_traverse_direct_children(gn, gtf_visitor,
                                                 save_exon_node, env);
  return has_err;
}

static int gtf_show_genome_feature(GenomeNode *gn, void *data, Env *env)
{
  GTFVisitor *gtf_visitor = (GTFVisitor*) data;
  GenomeFeatureType gft;
  int has_err = 0;
  switch ((gft = genome_feature_get_type((GenomeFeature*) gn))) {
    case gft_gene:
      gtf_visitor->gene_id++;
      gtf_visitor->transcript_id = 0;
      has_err = gtf_show_transcript(gn, gtf_visitor, env);
      break;
    case gft_mRNA:
      has_err = gtf_show_transcript(gn, gtf_visitor, env);
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
  return has_err;
}

static int gtf_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                      Env *env)
{
  GTFVisitor *gtf_visitor;
  int has_err;
  env_error_check(env);
  gtf_visitor = gtf_visitor_cast(gv);
  has_err = genome_node_traverse_children((GenomeNode*) gf, gtf_visitor,
                                          gtf_show_genome_feature, false, env);
  return has_err;
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

GenomeVisitor* gtf_visitor_new(GenFile *outfp, Env *env)
{
  GenomeVisitor *gv = genome_visitor_create(gtf_visitor_class(), env);
  GTFVisitor *gtf_visitor = gtf_visitor_cast(gv);
  gtf_visitor->gene_id = 0;
  gtf_visitor->exon_features = array_new(sizeof (GenomeNode*), env);
  gtf_visitor->outfp = outfp;
  return gv;
}
