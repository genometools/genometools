/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "libgtcore/queue.h"
#include "libgtcore/undef.h"
#include "libgtext/filter_visitor.h"
#include "libgtext/genome_visitor_rep.h"

struct FilterVisitor {
  const GenomeVisitor parent_instance;
  Queue *genome_node_buffer;
  Str *seqid,
      *typefilter;
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
  queue_delete(filter_visitor->genome_node_buffer);
  str_delete(filter_visitor->seqid);
  str_delete(filter_visitor->typefilter);
}

static int filter_visitor_comment(GenomeVisitor *gv, Comment *c, Error *e)
{
  FilterVisitor *filter_visitor;
  error_check(e);
  filter_visitor = filter_visitor_cast(gv);
  queue_add(filter_visitor->genome_node_buffer, c);
  return 0;
}

static int filter_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                         Error *e)
{
  FilterVisitor *fv;
  bool filter_node = false;
  error_check(e);
  fv = filter_visitor_cast(gv);
  if (!str_length(fv->seqid) || /* no seqid was specified or seqids are equal */
      !str_cmp(fv->seqid, genome_node_get_seqid((GenomeNode*) gf))) {
    /* enforce maximum gene length */
    /* XXX: we (spuriously) assume that genes are always root nodes */
    if (gf && genome_feature_get_type(gf) == gft_gene) {
      if (fv->max_gene_length != UNDEF_ULONG &&
          range_length(genome_node_get_range((GenomeNode*) gf)) >
          fv->max_gene_length) {
        filter_node = true;
      }
      else if (fv->max_gene_num != UNDEF_ULONG &&
               fv->gene_num >= fv->max_gene_num) {
        filter_node = true;
      }
      else if (fv->min_gene_score != UNDEF_DOUBLE &&
               genome_feature_get_score(gf) < fv->min_gene_score) {
        filter_node = true;
      }
      if (!filter_node)
        fv->gene_num++; /* gene passed filter */
    }
  }
  else
    filter_node = true;

  if (filter_node)
    genome_node_rec_delete((GenomeNode*) gf);
  else
    queue_add(fv->genome_node_buffer, gf);

  return 0;
}

static int filter_visitor_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                          Error *e)
{
  FilterVisitor *filter_visitor;
  error_check(e);
  filter_visitor = filter_visitor_cast(gv);
  if (!str_length(filter_visitor->seqid) || /* no seqid was specified */
      !str_cmp(filter_visitor->seqid,    /* or seqids are equal */
               genome_node_get_seqid((GenomeNode*) sr))) {
    queue_add(filter_visitor->genome_node_buffer, sr);
  }
  else
    genome_node_rec_delete((GenomeNode*) sr);
  return 0;
}

const GenomeVisitorClass* filter_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (FilterVisitor),
                                          filter_visitor_free,
                                          filter_visitor_comment,
                                          filter_visitor_genome_feature,
                                          filter_visitor_sequence_region };
  return &gvc;
}

GenomeVisitor* filter_visitor_new(Str *seqid, Str *typefilter,
                                  unsigned long max_gene_length,
                                  unsigned long max_gene_num,
                                  double min_gene_score)
{
  GenomeVisitor *gv = genome_visitor_create(filter_visitor_class());
  FilterVisitor *filter_visitor = filter_visitor_cast(gv);
  filter_visitor->genome_node_buffer = queue_new();
  filter_visitor->seqid = str_ref(seqid);
  filter_visitor->typefilter = str_ref(typefilter);
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
  FilterVisitor *filter_visitor;
  filter_visitor = filter_visitor_cast(gv);
  return queue_get(filter_visitor->genome_node_buffer);
}
