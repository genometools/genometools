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

#include "libgtext/genome_node.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtext/region_mapping.h"
#include "libgtext/reverse.h"

static int extract_join_feature(GenomeNode *gn, GenomeFeatureType type,
                                RegionMapping *region_mapping, Str *sequence,
                                bool *reverse_strand, Error *err)
{
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  GenomeFeature *gf;
  Range range;
  int had_err = 0;

  error_check(err);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);

  if (genome_feature_get_type(gf) == type) {
    had_err = region_mapping_get_raw_sequence(region_mapping, &raw_sequence,
                                              genome_node_get_seqid(gn), err);
    if (!had_err) {
      range = genome_node_get_range(gn);
      assert(range.start); /* 1-based coordinates */
      raw_sequence += range.start - 1;
      had_err = region_mapping_get_raw_sequence_length(region_mapping,
                                                       &raw_sequence_length,
                                                      genome_node_get_seqid(gn),
                                                       err);
    }
    if (!had_err) {
      assert(range.end <= raw_sequence_length);
      str_append_cstr_nt(sequence, raw_sequence, range_length(range));
      if (genome_feature_get_strand(gf) == STRAND_REVERSE)
        *reverse_strand = true;
    }
  }
  return had_err;
}

int extract_feat_sequence(Str *sequence, GenomeNode *gn, GenomeFeatureType type,
                          bool join, RegionMapping *region_mapping, Error *err)
{
  GenomeFeature *gf;
  Range range;
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  int had_err = 0;

  error_check(err);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);

  if (join) {
    GenomeNodeIterator *gni;
    GenomeNode *child;
    bool reverse_strand = false;
    /* in this case we have to traverse the children */
    gni = genome_node_iterator_new_direct(gn);
    while (!had_err && (child = genome_node_iterator_next(gni))) {
      if (extract_join_feature(child, type, region_mapping, sequence,
                               &reverse_strand, err)) {
        had_err = -1;
      }
    }
    genome_node_iterator_delete(gni);
    if (!had_err && str_length(sequence)) {
      if (reverse_strand) {
        had_err = reverse_complement(str_get(sequence),
                                     str_length(sequence), err);
      }
    }
  }
  else if (genome_feature_get_type(gf) == type) {
    assert(!had_err);
    /* otherwise we only have to look this feature */
    range = genome_node_get_range(gn);
    assert(range.start); /* 1-based coordinates */
    had_err = region_mapping_get_raw_sequence_length(region_mapping,
                                                     &raw_sequence_length,
                                                     genome_node_get_seqid(gn),
                                                     err);
    if (!had_err) {
      assert(range.end <= raw_sequence_length);
      had_err = region_mapping_get_raw_sequence(region_mapping,
                                                &raw_sequence,
                                                genome_node_get_seqid(gn), err);
    }
    if (!had_err) {
      str_append_cstr_nt(sequence, raw_sequence + range.start - 1,
                         range_length(range));
      if (genome_feature_get_strand(gf) == STRAND_REVERSE) {
        had_err = reverse_complement(str_get(sequence), str_length(sequence),
                                     err);
      }
    }
  }
  return had_err;
}
