/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "bioseq.h"
#include "error.h"
#include "extractfeat_visitor.h"
#include "fasta.h"
#include "genome_visitor_rep.h"
#include "reverse.h"
#include "translate.h"

struct ExtractFeatVisitor {
  const GenomeVisitor parent_instance;
  Str *sequence_file, /* the (current) sequence file */
      *description,   /* the description of the currently extracted feature */
      *sequence,      /* the sequence of the currently extracted feature */
      *protein;       /* the translated protein sequence (if applicable) */
  GenomeFeatureType type;
  bool join,
       translate,
       reverse_strand;
  Bioseq *bioseq;
  unsigned long fastaseq_counter;
  RegionMapping *regionmapping;
};

#define extractfeat_visitor_cast(GV)\
        genome_visitor_cast(extractfeat_visitor_class(), GV)

static void extractfeat_visitor_free(GenomeVisitor *gv)
{
  ExtractFeatVisitor *extractfeat_visitor = extractfeat_visitor_cast(gv);
  assert(extractfeat_visitor);
  str_free(extractfeat_visitor->sequence_file);
  str_free(extractfeat_visitor->description);
  str_free(extractfeat_visitor->sequence);
  str_free(extractfeat_visitor->protein);
  bioseq_free(extractfeat_visitor->bioseq);
  regionmapping_free(extractfeat_visitor->regionmapping);
}

static int extract_join_feature(GenomeNode *gn, void *data, Error *err)
{
  ExtractFeatVisitor *v = (ExtractFeatVisitor*) data;
  GenomeFeatureType gf_type;
  const char *raw_sequence;
  Genome_feature *gf;
  Range range;

  error_check(err);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  gf_type = genome_feature_get_type(gf);

  if (gf_type == v->type) {
    raw_sequence = bioseq_get_raw_sequence(v->bioseq);
    range = genome_node_get_range(gn);
    assert(range.start); /* 1-based coordinates */
    raw_sequence += range.start - 1;
    assert(range.end <= bioseq_get_raw_sequence_length(v->bioseq));
    str_append_cstr_nt(v->sequence, raw_sequence, range_length(range));
    if (genome_feature_get_strand(gf) == STRAND_REVERSE)
      v->reverse_strand = true;
  }
  return 0;
}

static int extract_feature(GenomeNode *gn, void *data, Error *err)
{
  ExtractFeatVisitor *v = (ExtractFeatVisitor*) data;
  GenomeFeatureType gf_type;
  Genome_feature *gf;
  Range range;
  int has_err = 0;

  error_check(err);
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
    has_err = genome_node_traverse_direct_children(gn, v, extract_join_feature,
                                                   err);
    if (!has_err && str_length(v->sequence)) {
      if (v->reverse_strand) {
        has_err = reverse_complement(str_get(v->sequence),
                                     str_length(v->sequence), err);
      }
      if (!has_err) {
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
    assert(!has_err);
    /* otherwise we only have to look this feature */
    range = genome_node_get_range(gn);
    assert(range.start); /* 1-based coordinates */
    assert(range.end <= bioseq_get_raw_sequence_length(v->bioseq));
    str_reset(v->sequence);
    str_append_cstr_nt(v->sequence, bioseq_get_raw_sequence(v->bioseq) +
                       range.start - 1, range_length(range));
    if (genome_feature_get_strand(gf) == STRAND_REVERSE) {
      has_err = reverse_complement(str_get(v->sequence),
                                   str_length(v->sequence), err);
    }
    if (!has_err) {
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
  return has_err;
}

static int extractfeat_visitor_genome_feature(GenomeVisitor *gv,
                                               Genome_feature *gf,
                                               /*@unused@*/ Log *l, Error *err)
{
  ExtractFeatVisitor *efv;
  error_check(err);
  efv = extractfeat_visitor_cast(gv);
  if (efv->regionmapping) { /* region mapping used -> determine bioseq */
    if (!efv->sequence_file ||
        str_cmp(efv->sequence_file, genome_node_get_seqid((GenomeNode*) gf))) {
      str_free(efv->sequence_file);
      efv->sequence_file = regionmapping_map(efv->regionmapping,
                              str_get(genome_node_get_seqid((GenomeNode*) gf)));
      bioseq_free(efv->bioseq);
      efv->bioseq = bioseq_new_str(efv->sequence_file);
    }
  }
  return genome_node_traverse_children((GenomeNode*) gf, efv, extract_feature,
                                       false, err);
}

static int extractfeat_visitor_sequence_region(GenomeVisitor *gv,
                                               SequenceRegion *sr,
                                               /*@unused@*/ Log *l, Error *err)
{
  ExtractFeatVisitor *efv;
  int has_err = 0;
  error_check(err);
  efv = extractfeat_visitor_cast(gv);
  if (!efv->regionmapping) { /* we have only one bioseq */
    /* check if the given sequence file contains this sequence (region) */
    if (!bioseq_contains_sequence(efv->bioseq,
                                  str_get(genome_node_get_seqid((GenomeNode*)
                                                              sr)))) {
      error_set(err, "sequence \"%s\" not contained in sequence file \"%s\"",
                str_get(genome_node_get_seqid((GenomeNode*) sr)),
                str_get(efv->sequence_file));
      has_err = -1;
    }
  }
  return has_err;
}

const GenomeVisitorClass* extractfeat_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof(ExtractFeatVisitor),
                                          extractfeat_visitor_free,
                                          NULL,
                                          extractfeat_visitor_genome_feature,
                                          extractfeat_visitor_sequence_region,
                                          NULL };
  return &gvc;
}

static GenomeVisitor* extractfeat_visitor_new(GenomeFeatureType type,
                                              bool join, bool translate)
{
  GenomeVisitor *gv = genome_visitor_create(extractfeat_visitor_class());
  ExtractFeatVisitor *efv= extractfeat_visitor_cast(gv);
  efv->description = str_new();
  efv->sequence = str_new();
  efv->protein = str_new();
  efv->type = type;
  efv->join = join;
  efv->translate = translate;
  efv->fastaseq_counter = 0;
  return gv;
}

GenomeVisitor* extractfeat_visitor_new_seqfile(Str *sequence_file,
                                               GenomeFeatureType type,
                                               bool join, bool translate)
{
  GenomeVisitor *gv;
  ExtractFeatVisitor *efv;
  assert(sequence_file);
  gv = extractfeat_visitor_new(type, join, translate);
  efv = extractfeat_visitor_cast(gv);
  efv->sequence_file = str_ref(sequence_file);
  efv->bioseq = bioseq_new_str(sequence_file);
  return gv;
}

GenomeVisitor* extractfeat_visitor_new_regionmapping(RegionMapping *rm,
                                                     GenomeFeatureType type,
                                                     bool join, bool translate)
{
  GenomeVisitor *gv;
  ExtractFeatVisitor *efv;
  gv = extractfeat_visitor_new(type, join, translate);
  efv = extractfeat_visitor_cast(gv);
  efv->regionmapping = rm;
  return gv;
}
