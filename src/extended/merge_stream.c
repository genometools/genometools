/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/ma.h"
#include "core/undef.h"
#include "extended/genome_stream_rep.h"
#include "extended/merge_stream.h"
#include "extended/sequence_region.h"

struct MergeStream {
  const GenomeStream parent_instance;
  GT_Array *genome_streams;
  GT_GenomeNode **buffer;
};

#define merge_stream_cast(GS)\
        genome_stream_cast(merge_stream_class(), GS)

int merge_stream_next_tree(GenomeStream *gs, GT_GenomeNode **gn, GT_Error *err)
{
  MergeStream *ms;
  GT_GenomeNode *min_node = NULL;
  unsigned long i, j, min_i = UNDEF_ULONG;
  unsigned int genome_node_consolidated;
  int had_err = 0;

  gt_error_check(err);

  ms = merge_stream_cast(gs);

  /* fill buffers */
  for (i = 0; i < gt_array_size(ms->genome_streams); i++) {
    if (!ms->buffer[i]) {
      had_err = genome_stream_next_tree(*(GenomeStream**)
                                        gt_array_get(ms->genome_streams, i),
                                        ms->buffer + i, err);
      if (had_err)
        break;
    }
  }

  /* consolidate sequence regions (to avoid duplicates) */
  if (!had_err) {
    for (;;) {
      genome_node_consolidated = 0;
      for (i = 0; i < gt_array_size(ms->genome_streams); i++) {
        for (j = i+1; j < gt_array_size(ms->genome_streams); j++) {
          assert(i != j);
          if (genome_nodes_are_equal_sequence_regions(ms->buffer[i],
                                                      ms->buffer[j])) {
            sequence_regions_consolidate(ms->buffer[i], ms->buffer[j]);
            genome_node_rec_delete(ms->buffer[j]);
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
    for (i = 0; i < gt_array_size(ms->genome_streams); i++) {
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

static void merge_stream_free(GenomeStream *gs)
{
  MergeStream *ms = merge_stream_cast(gs);
  unsigned long i;
  for (i = 0; i < gt_array_size(ms->genome_streams); i++)
    genome_stream_delete(*(GenomeStream**) gt_array_get(ms->genome_streams, i));
  gt_array_delete(ms->genome_streams);
  ma_free(ms->buffer);
}

const GenomeStreamClass* merge_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (MergeStream),
                                         merge_stream_next_tree,
                                         merge_stream_free };
  return &gsc;
}

GenomeStream* merge_stream_new(const GT_Array *genome_streams)
{
  GenomeStream *in_stream,
               *gs = genome_stream_create(merge_stream_class(), true);
  MergeStream *ms = merge_stream_cast(gs);
  unsigned long i;
#ifndef NDEBUG
  assert(gt_array_size(genome_streams)); /* at least on input stream given */
  /* each input stream is sorted */
  for (i = 0; i < gt_array_size(genome_streams); i++) {
    assert(genome_stream_is_sorted(*(GenomeStream**)
                                   gt_array_get(genome_streams, i)));
  }
#endif
  ms->genome_streams = gt_array_new(sizeof (GenomeStream*));
  for (i = 0; i < gt_array_size(genome_streams); i++) {
    in_stream = genome_stream_ref(*(GenomeStream**)
                                  gt_array_get(genome_streams, i));
    gt_array_add(ms->genome_streams, in_stream);
  }
  ms->buffer = ma_calloc(gt_array_size(genome_streams), sizeof (GT_GenomeNode*));
  return gs;
}
