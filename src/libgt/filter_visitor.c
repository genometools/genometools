/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "filter_visitor.h"
#include "genome_visitor_rep.h"
#include "queue.h"
#include "undef.h"

struct FilterVisitor {
  const GenomeVisitor parent_instance;
  Queue *genome_node_buffer;
  Str *seqid;
  unsigned long max_gene_length,
                gene_num,     /* the number of passed genes */
                max_gene_num; /* the maximal number of genes which can pass */
  double min_gene_score;
};

#define filter_visitor_cast(GV)\
        genome_visitor_cast(filter_visitor_class(), GV)

static void filter_visitor_free(GenomeVisitor *gv)
{
  FilterVisitor *filter_visitor = filter_visitor_cast(gv);
  queue_free(filter_visitor->genome_node_buffer);
  str_free(filter_visitor->seqid);
}

static int filter_visitor_comment(GenomeVisitor *gv, Comment *c,
                                  Log *l, Error *err)
{
  FilterVisitor *filter_visitor;
  error_check(err);
  filter_visitor = filter_visitor_cast(gv);
  queue_add(filter_visitor->genome_node_buffer, c);
  return 0;
}

static int filter_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                         Log *l, Error *err)
{
  FilterVisitor *filter_visitor;
  error_check(err);
  filter_visitor = filter_visitor_cast(gv);
  if (!str_get(filter_visitor->seqid) || /* no seqid was specified */
      !str_cmp(filter_visitor->seqid,    /* or seqids are equal */
               genome_node_get_seqid((GenomeNode*) gf))) {
    /* enforce maximum gene length */
    /* XXX: we (spuriously) assume that genes are always root nodes */
    if (gf && genome_feature_get_type(gf) == gft_gene) {
      if (filter_visitor->max_gene_length != UNDEFULONG &&
          range_length(genome_node_get_range((GenomeNode*) gf)) >
          filter_visitor->max_gene_length) {
        return 0;
      }
      else if (filter_visitor->max_gene_num != UNDEFULONG &&
               filter_visitor->gene_num >= filter_visitor->max_gene_num) {
        return 0;
      }
      else if (filter_visitor->min_gene_score != UNDEFDOUBLE &&
               genome_feature_get_score(gf) < filter_visitor->min_gene_score) {
        return 0;
      }
      filter_visitor->gene_num++; /* gene passed filter */
    }
    queue_add(filter_visitor->genome_node_buffer, gf);
  }
  return 0;
}

static int filter_visitor_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                          Log *l, Error *err)
{
  FilterVisitor *filter_visitor;
  error_check(err);
  filter_visitor = filter_visitor_cast(gv);
  if (!str_get(filter_visitor->seqid) || /* no seqid was specified */
      !str_cmp(filter_visitor->seqid,    /* or seqids are equal */
               genome_node_get_seqid((GenomeNode*) sr))) {
    queue_add(filter_visitor->genome_node_buffer, sr);
  }
  return 0;
}

const GenomeVisitorClass* filter_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof(FilterVisitor),
                                          filter_visitor_free,
                                          filter_visitor_comment,
                                          filter_visitor_genome_feature,
                                          filter_visitor_sequence_region,
                                          NULL };
  return &gvc;
}

GenomeVisitor* filter_visitor_new(Str *seqid, unsigned long max_gene_length,
                                  unsigned long max_gene_num,
                                  double min_gene_score)
{
  GenomeVisitor *gv = genome_visitor_create(filter_visitor_class());
  FilterVisitor *filter_visitor = filter_visitor_cast(gv);
  filter_visitor->genome_node_buffer = queue_new(sizeof(GenomeNode*));
  filter_visitor->seqid = str_ref(seqid);
  filter_visitor->max_gene_length = max_gene_length;
  filter_visitor->gene_num = 0;
  filter_visitor->max_gene_num = max_gene_num;
  filter_visitor->min_gene_score = min_gene_score;
  return gv;
}

unsigned long filter_visitor_node_buffer_size(GenomeVisitor *gv)
{
  FilterVisitor *filter_visitor = filter_visitor_cast(gv);
  return queue_size(filter_visitor->genome_node_buffer);
}

GenomeNode* filter_visitor_get_node(GenomeVisitor *gv)
{
  FilterVisitor *filter_visitor = filter_visitor_cast(gv);
  return *(GenomeNode**) queue_get(filter_visitor->genome_node_buffer);
}
