/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "libgtcore/undef.h"
#include <libgtext/genome_stream_rep.h>
#include <libgtext/merge_stream.h>

struct MergeStream {
  const GenomeStream parent_instance;
  Array *genome_streams;
  GenomeNode **buffer;
};

#define merge_stream_cast(GS)\
        genome_stream_cast(merge_stream_class(), GS)

static unsigned int genome_nodes_are_equal_sequence_regions(GenomeNode *gn_a,
                                                            GenomeNode *gn_b)
{
  void *sr_a, *sr_b;

  sr_a = gn_a ? genome_node_cast(sequence_region_class(), gn_a) : NULL;
  sr_b = gn_b ? genome_node_cast(sequence_region_class(), gn_b) : NULL;

  if (sr_a && sr_b &&
      str_cmp(genome_node_get_seqid(gn_a), genome_node_get_seqid(gn_b)) == 0) {
    return 1;
  }
  return 0;
}

static void consolidate_sequence_regions(GenomeNode *gn_a, GenomeNode *gn_b)
{
  Range range_a, range_b;

  assert(gn_a);
  assert(gn_b);
  assert(genome_node_cast(sequence_region_class(), gn_a));
  assert(genome_node_cast(sequence_region_class(), gn_b));
  assert(str_cmp(genome_node_get_seqid(gn_a),
                 genome_node_get_seqid(gn_b)) == 0);

  range_a = genome_node_get_range(gn_a);
  range_b = genome_node_get_range(gn_b);
  range_a = range_join(range_a, range_b);
  genome_node_set_range(gn_a, range_a);
}

int merge_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  MergeStream *ms;
  GenomeNode *min_node = NULL;
  unsigned long i, j, min_i = UNDEF_ULONG;
  unsigned int genome_node_consolidated;
  int had_err = 0;

  env_error_check(env);

  ms = merge_stream_cast(gs);

  /* fill buffers */
  for (i = 0; i < array_size(ms->genome_streams); i++) {
    if (!ms->buffer[i]) {
      had_err = genome_stream_next_tree(*(GenomeStream**)
                                        array_get(ms->genome_streams, i),
                                        ms->buffer + i, env);
      if (had_err)
        break;
    }
  }

  /* consolidate sequence regions (to avoid duplicates) */
  if (!had_err) {
    for (;;) {
      genome_node_consolidated = 0;
      for (i = 0; i < array_size(ms->genome_streams); i++) {
        for (j = i+1; j < array_size(ms->genome_streams); j++) {
          assert(i != j);
          if (genome_nodes_are_equal_sequence_regions(ms->buffer[i],
                                                      ms->buffer[j])) {
            consolidate_sequence_regions(ms->buffer[i], ms->buffer[j]);
            genome_node_rec_delete(ms->buffer[j], env);
            ms->buffer[j] = NULL;
          }
        }
      }
      if (!genome_node_consolidated)
        break;
    }
  }

  /* find minimal node */
  if (!had_err) {
    for (i = 0; i < array_size(ms->genome_streams); i++) {
      if (ms->buffer[i]) {
        if (min_i != UNDEF_ULONG) {
          if (genome_node_compare(ms->buffer + i, ms->buffer + min_i) < 0)
            min_i = i;
        }
        else min_i = i;
      }
    }
    if (min_i != UNDEF_ULONG) {
      min_node = ms->buffer[min_i];
      ms->buffer[min_i] = NULL;
    }
  }

  *gn = min_node;
  return had_err;
}

static void merge_stream_free(GenomeStream *gs, Env *env)
{
  MergeStream *ms = merge_stream_cast(gs);
  array_delete(ms->genome_streams, env);
  env_ma_free(ms->buffer, env);
}

const GenomeStreamClass* merge_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (MergeStream),
                                         merge_stream_next_tree,
                                         merge_stream_free };
  return &gsc;
}

GenomeStream* merge_stream_new(const Array *genome_streams, Env *env)
{
  GenomeStream *gs = genome_stream_create(merge_stream_class(), true, env);
  MergeStream *ms = merge_stream_cast(gs);
#ifndef NDEBUG
  unsigned long i;
  assert(array_size(genome_streams)); /* at least on input stream given */
  /* each input stream is sorted */
  for (i = 0; i < array_size(genome_streams); i++) {
    assert(genome_stream_is_sorted(*(GenomeStream**)
                                   array_get(genome_streams, i)));
  }
#endif
  ms->genome_streams = array_clone(genome_streams, env);
  ms->buffer = env_ma_calloc(env, array_size(genome_streams),
                             sizeof (GenomeNode*));
  return gs;
}
