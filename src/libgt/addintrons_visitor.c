/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "addintrons_visitor.h"
#include "genome_visitor_rep.h"
#include "undef.h"

struct AddIntronsVisitor {
  const GenomeVisitor parent_instance;
  GenomeFeature *parent_feature,
                *previous_exon_feature;
};

#define addintrons_visitor_cast(GV)\
        genome_visitor_cast(addintrons_visitor_class(), GV)

static void addintrons_visitor_free(GenomeVisitor *gv, Env *env)
{
  AddIntronsVisitor *addintrons_visitor = addintrons_visitor_cast(gv);
  assert(addintrons_visitor);
}

static int addintrons_in_children(GenomeNode *gn, void *data, Env *env)
{
  AddIntronsVisitor *v = (AddIntronsVisitor*) data;
  GenomeFeature *current_feature;
  GenomeNode *intron_node;
  Range previous_range, current_range, intron_range;
  Strand previous_strand, current_strand, intron_strand;
  Str *parent_seqid;
  env_error_check(env);
  current_feature = genome_node_cast(genome_feature_class(), gn);
  assert(current_feature);
  if (genome_feature_get_type(current_feature) == gft_exon) {
    if (v->previous_exon_feature) {
      /* determine intron range */
      previous_range = genome_node_get_range((GenomeNode*)
                                             v->previous_exon_feature);
      current_range = genome_node_get_range(gn);
      assert(previous_range.end < current_range.start);
      intron_range.start = previous_range.end + 1;
      intron_range.end = current_range.start - 1;

      /* determine intron strand */
      previous_strand = genome_feature_get_strand(v->previous_exon_feature);
      current_strand = genome_feature_get_strand(current_feature);
      assert(previous_strand == current_strand);
      intron_strand = previous_strand;

      /* determine sequence id */
      parent_seqid = genome_node_get_seqid((GenomeNode*) v->parent_feature);
      assert(!str_cmp(parent_seqid,
             genome_node_get_seqid((GenomeNode*) v->previous_exon_feature)));
      assert(!str_cmp(parent_seqid,
             genome_node_get_seqid((GenomeNode*) current_feature)));

      /* create intron */
      intron_node = genome_feature_new(gft_intron, intron_range, intron_strand,
                                       NULL, UNDEFULONG, env);
      genome_node_set_seqid(intron_node, parent_seqid);
      genome_node_is_part_of_genome_node((GenomeNode*) v->parent_feature,
                                         intron_node, env);
    }
    v->previous_exon_feature = current_feature;
  }
  return 0;
}

static int addintrons_if_necessary(GenomeNode *gn, void *data, Env *env)
{
  AddIntronsVisitor *v = (AddIntronsVisitor*) data;
  GenomeFeature *gf;
  env_error_check(env);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  v->parent_feature = gf;
  v->previous_exon_feature = NULL;
  return genome_node_traverse_direct_children(gn, v, addintrons_in_children,
                                              env);
}

static int addintrons_visitor_genome_feature(GenomeVisitor *gv,
                                             GenomeFeature *gf, Env *env)
{
  AddIntronsVisitor *v;
  env_error_check(env);
  v = addintrons_visitor_cast(gv);
  return genome_node_traverse_children((GenomeNode*) gf, v,
                                       addintrons_if_necessary, false, env);
}

const GenomeVisitorClass* addintrons_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (AddIntronsVisitor),
                                          addintrons_visitor_free,
                                          NULL,
                                          addintrons_visitor_genome_feature,
                                          NULL,
                                          NULL };
  return &gvc;
}

GenomeVisitor* addintrons_visitor_new(Env *env)
{
  return genome_visitor_create(addintrons_visitor_class(), env);
}
