/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/extractfeat_visitor.h"
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
  RegionMapping *regionmapping;
};

#define extractfeat_visitor_cast(GV)\
        genome_visitor_cast(extractfeat_visitor_class(), GV)

static void extractfeat_visitor_free(GenomeVisitor *gv)
{
  ExtractFeatVisitor *extractfeat_visitor = extractfeat_visitor_cast(gv);
  assert(extractfeat_visitor);
  str_delete(extractfeat_visitor->description);
  str_delete(extractfeat_visitor->sequence);
  str_delete(extractfeat_visitor->protein);
  regionmapping_delete(extractfeat_visitor->regionmapping);
}

static int extract_join_feature(GenomeNode *gn, void *data, Error *e)
{
  ExtractFeatVisitor *v = (ExtractFeatVisitor*) data;
  GenomeFeatureType gf_type;
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  GenomeFeature *gf;
  Range range;
  int had_err = 0;

  error_check(e);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  gf_type = genome_feature_get_type(gf);

  if (gf_type == v->type) {
    had_err = regionmapping_get_raw_sequence(v->regionmapping, &raw_sequence,
                                             genome_node_get_seqid(gn), e);
    if (!had_err) {
      range = genome_node_get_range(gn);
      assert(range.start); /* 1-based coordinates */
      raw_sequence += range.start - 1;
      had_err = regionmapping_get_raw_sequence_length(v->regionmapping,
                                                      &raw_sequence_length,
                                                      genome_node_get_seqid(gn),
                                                      e);
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

static int extract_feature(GenomeNode *gn, void *data, Error *e)
{
  ExtractFeatVisitor *v = (ExtractFeatVisitor*) data;
  GenomeFeatureType gf_type;
  GenomeFeature *gf;
  Range range;
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  int had_err = 0;

  error_check(e);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  gf_type = genome_feature_get_type(gf);

  /* construct description if necessary */
  if (!str_length(v->description)) {
    str_append_cstr(v->description, genome_feature_type_get_cstr(v->type));
    str_append_char(v->description, '_');
    v->fastaseq_counter++;
    str_append_ulong(v->description, v->fastaseq_counter);
    if (v->join)
      str_append_cstr(v->description, " (joined)");
    if (v->translate)
      str_append_cstr(v->description, " (translated)");
  }

  if (v->join) {
    /* in this case we have to traverse the children */
    str_reset(v->sequence);
    v->reverse_strand = false;
    had_err = genome_node_traverse_direct_children(gn, v, extract_join_feature,
                                                   e);
    if (!had_err && str_length(v->sequence)) {
      if (v->reverse_strand) {
        had_err = reverse_complement(str_get(v->sequence),
                                     str_length(v->sequence), e);
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
    had_err = regionmapping_get_raw_sequence_length(v->regionmapping,
                                                    &raw_sequence_length,
                                                    genome_node_get_seqid(gn),
                                                    e);
    if (!had_err) {
      assert(range.end <= raw_sequence_length);
      str_reset(v->sequence);
      had_err = regionmapping_get_raw_sequence(v->regionmapping, &raw_sequence,
                                               genome_node_get_seqid(gn), e);
    }
    if (!had_err) {
      str_append_cstr_nt(v->sequence, raw_sequence + range.start - 1,
                         range_length(range));
      if (genome_feature_get_strand(gf) == STRAND_REVERSE) {
        had_err = reverse_complement(str_get(v->sequence),
                                     str_length(v->sequence), e);
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

static int extractfeat_visitor_genome_feature(GenomeVisitor *gv,
                                              GenomeFeature *gf, Error *e)
{
  ExtractFeatVisitor *efv;
  error_check(e);
  efv = extractfeat_visitor_cast(gv);
  assert(efv->regionmapping);
  return genome_node_traverse_children((GenomeNode*) gf, efv, extract_feature,
                                       false, e);
}

const GenomeVisitorClass* extractfeat_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (ExtractFeatVisitor),
                                          extractfeat_visitor_free,
                                          NULL,
                                          extractfeat_visitor_genome_feature,
                                          NULL };
  return &gvc;
}

GenomeVisitor* extractfeat_visitor_new(RegionMapping *rm,
                                       GenomeFeatureType type, bool join,
                                       bool translate)
{
  GenomeVisitor *gv;
  ExtractFeatVisitor *efv;
  assert(rm);
  gv = genome_visitor_create(extractfeat_visitor_class());
  efv= extractfeat_visitor_cast(gv);
  efv->description = str_new();
  efv->sequence = str_new();
  efv->protein = str_new();
  efv->type = type;
  efv->join = join;
  efv->translate = translate;
  efv->fastaseq_counter = 0;
  efv->regionmapping = rm;
  return gv;
}
