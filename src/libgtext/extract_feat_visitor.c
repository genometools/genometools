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
#include "libgtcore/fasta.h"
#include "libgtcore/translate.h"
#include "libgtext/extract_feat_visitor.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtext/genome_visitor_rep.h"
#include "libgtext/reverse.h"

struct ExtractFeatVisitor {
  const GenomeVisitor parent_instance;
  Str *description,   /* the description of the currently extracted feature */
      *sequence,      /* the sequence of the currently extracted feature */
      *protein;       /* the translated protein sequence (if applicable) */
  GenomeFeatureType type;
  bool join,
       translate,
       reverse_strand;
  unsigned long fastaseq_counter;
  RegionMapping *region_mapping;
};

#define extract_feat_visitor_cast(GV)\
        genome_visitor_cast(extract_feat_visitor_class(), GV)

static void extract_feat_visitor_free(GenomeVisitor *gv)
{
  ExtractFeatVisitor *extract_feat_visitor = extract_feat_visitor_cast(gv);
  assert(extract_feat_visitor);
  str_delete(extract_feat_visitor->description);
  str_delete(extract_feat_visitor->sequence);
  str_delete(extract_feat_visitor->protein);
  region_mapping_delete(extract_feat_visitor->region_mapping);
}

static int extract_join_feature(GenomeNode *gn, ExtractFeatVisitor *v,
                                Error *err)
{
  GenomeFeatureType gf_type;
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  GenomeFeature *gf;
  Range range;
  int had_err = 0;

  error_check(err);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  gf_type = genome_feature_get_type(gf);

  if (gf_type == v->type) {
    had_err = region_mapping_get_raw_sequence(v->region_mapping, &raw_sequence,
                                              genome_node_get_seqid(gn), err);
    if (!had_err) {
      range = genome_node_get_range(gn);
      assert(range.start); /* 1-based coordinates */
      raw_sequence += range.start - 1;
      had_err = region_mapping_get_raw_sequence_length(v->region_mapping,
                                                       &raw_sequence_length,
                                                      genome_node_get_seqid(gn),
                                                       err);
    }
    if (!had_err) {
      assert(range.end <= raw_sequence_length);
      str_append_cstr_nt(v->sequence, raw_sequence, range_length(range));
      if (genome_feature_get_strand(gf) == STRAND_REVERSE)
        v->reverse_strand = true;
    }
  }
  return had_err;
}

static void construct_description(Str *description, GenomeFeatureType type,
                                  unsigned long counter, bool join,
                                  bool translate)
{
  assert(!str_length(description));
  str_append_cstr(description, genome_feature_type_get_cstr(type));
  str_append_char(description, '_');
  str_append_ulong(description, counter);
  if (join)
    str_append_cstr(description, " (joined)");
  if (translate)
    str_append_cstr(description, " (translated)");
}

static int extract_feature(GenomeNode *gn, ExtractFeatVisitor *v, Error *err)
{
  GenomeFeatureType gf_type;
  GenomeFeature *gf;
  Range range;
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  int had_err = 0;

  error_check(err);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  gf_type = genome_feature_get_type(gf);

  /* construct description if necessary */
  if (!str_length(v->description)) {
    v->fastaseq_counter++;
    construct_description(v->description, v->type, v->fastaseq_counter, v->join,
                          v->translate);
  }

  if (v->join) {
    GenomeNodeIterator *gni;
    GenomeNode *child;
    /* in this case we have to traverse the children */
    str_reset(v->sequence);
    v->reverse_strand = false;
    gni = genome_node_iterator_new_direct(gn);
    while (!had_err && (child = genome_node_iterator_next(gni))) {
      if (extract_join_feature(child, v, err))
        had_err = -1;
    }
    genome_node_iterator_delete(gni);
    if (!had_err && str_length(v->sequence)) {
      if (v->reverse_strand) {
        had_err = reverse_complement(str_get(v->sequence),
                                     str_length(v->sequence), err);
      }
      if (!had_err) {
        if (v->translate) {
          str_reset(v->protein);
          translate_dna(v->protein, str_get(v->sequence),
                        str_length(v->sequence), 0);
          fasta_show_entry(str_get(v->description), str_get(v->protein),
                           str_length(v->protein), 0);
        }
        else {
          fasta_show_entry(str_get(v->description), str_get(v->sequence),
                           str_length(v->sequence), 0);
        }
      }
      str_reset(v->description);
    }
  }
  else if (gf_type == v->type) {
    assert(!had_err);
    /* otherwise we only have to look this feature */
    range = genome_node_get_range(gn);
    assert(range.start); /* 1-based coordinates */
    had_err = region_mapping_get_raw_sequence_length(v->region_mapping,
                                                     &raw_sequence_length,
                                                     genome_node_get_seqid(gn),
                                                     err);
    if (!had_err) {
      assert(range.end <= raw_sequence_length);
      str_reset(v->sequence);
      had_err = region_mapping_get_raw_sequence(v->region_mapping,
                                                &raw_sequence,
                                                genome_node_get_seqid(gn), err);
    }
    if (!had_err) {
      str_append_cstr_nt(v->sequence, raw_sequence + range.start - 1,
                         range_length(range));
      if (genome_feature_get_strand(gf) == STRAND_REVERSE) {
        had_err = reverse_complement(str_get(v->sequence),
                                     str_length(v->sequence), err);
      }
    }
    if (!had_err) {
      if (v->translate) {
        str_reset(v->protein);
        translate_dna(v->protein, str_get(v->sequence), str_length(v->sequence),
                      0);
        fasta_show_entry(str_get(v->description), str_get(v->protein),
                         str_length(v->protein), 0);
      }
      else {
        fasta_show_entry(str_get(v->description), str_get(v->sequence),
                         str_length(v->sequence), 0);
      }
    }
    str_reset(v->description);
  }
  return had_err;
}

static int extract_feat_visitor_genome_feature(GenomeVisitor *gv,
                                               GenomeFeature *gf, Error *err)
{
  ExtractFeatVisitor *efv;
  GenomeNodeIterator *gni;
  GenomeNode *gn;
  int had_err = 0;
  error_check(err);
  efv = extract_feat_visitor_cast(gv);
  assert(efv->region_mapping);
  gni = genome_node_iterator_new((GenomeNode*) gf);
  while (!had_err && (gn = genome_node_iterator_next(gni))) {
    if (extract_feature(gn, efv, err))
      had_err = -1;
  }
  genome_node_iterator_delete(gni);
  return had_err;
}

const GenomeVisitorClass* extract_feat_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (ExtractFeatVisitor),
                                          extract_feat_visitor_free,
                                          NULL,
                                          extract_feat_visitor_genome_feature,
                                          NULL };
  return &gvc;
}

GenomeVisitor* extract_feat_visitor_new(RegionMapping *rm,
                                       GenomeFeatureType type, bool join,
                                       bool translate)
{
  GenomeVisitor *gv;
  ExtractFeatVisitor *efv;
  assert(rm);
  gv = genome_visitor_create(extract_feat_visitor_class());
  efv= extract_feat_visitor_cast(gv);
  efv->description = str_new();
  efv->sequence = str_new();
  efv->protein = str_new();
  efv->type = type;
  efv->join = join;
  efv->translate = translate;
  efv->fastaseq_counter = 0;
  efv->region_mapping = rm;
  return gv;
}
