/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "filter_stream.h"
#include "genome_feature.h"
#include "genome_stream_rep.h"
#include "undef.h"

struct Filter_stream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  unsigned long max_gene_length;
  double min_gene_score;
};

#define filter_stream_cast(GS)\
        genome_stream_cast(filter_stream_class(), GS);

static int filter_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Log *l,
                                   Error *err)
{
  Filter_stream *fs;
  Genome_feature *gf;
  int has_err;

  error_check(err);

  fs = filter_stream_cast(gs);

  /* enforce maximum gene length */
  /* XXX: we (spuriously) assume that genes are always root nodes */
  while (!(has_err = genome_stream_next_tree(fs->in_stream, gn, l, err)) &&
         *gn) {
    /* XXX: use visitor pattern */
    gf = genome_node_cast(genome_feature_class(), *gn);
    if (gf && genome_feature_get_type(gf) == gft_gene) {
      if (fs->max_gene_length != UNDEFULONG &&
         range_length(genome_node_get_range(*gn)) > fs->max_gene_length) {
        genome_node_rec_free(*gn);
      }
      else if (fs->min_gene_score != UNDEFDOUBLE &&
               genome_feature_get_score(gf) < fs->min_gene_score) {
        genome_node_rec_free(*gn);
      }
      else
        break;
    }
    else
      break;
  }

  return has_err;
}

const GenomeStreamClass* filter_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (Filter_stream),
                                           filter_stream_next_tree,
                                           NULL };
  return &gsc;
}

GenomeStream* filter_stream_new(GenomeStream *in_stream,
                                 unsigned long max_gene_length,
                                 double min_gene_score)
{
  GenomeStream *gs = genome_stream_create(filter_stream_class(),
                                           genome_stream_is_sorted(in_stream));
  Filter_stream *filter_stream = filter_stream_cast(gs);
  assert(in_stream);
  filter_stream->in_stream = in_stream;
  filter_stream->max_gene_length = max_gene_length;
  filter_stream->min_gene_score = min_gene_score;
  return gs;
}
