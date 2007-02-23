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
      *sequence_name, /* the (current) sequence name */
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

static void extractfeat_visitor_free(GenomeVisitor *gv, Env *env)
{
  ExtractFeatVisitor *extractfeat_visitor = extractfeat_visitor_cast(gv);
  assert(extractfeat_visitor);
  str_delete(extractfeat_visitor->sequence_file, env);
  str_delete(extractfeat_visitor->sequence_name, env);
  str_delete(extractfeat_visitor->description, env);
  str_delete(extractfeat_visitor->sequence, env);
  str_delete(extractfeat_visitor->protein, env);
  bioseq_delete(extractfeat_visitor->bioseq, env);
  regionmapping_delete(extractfeat_visitor->regionmapping, env);
}

static int extract_join_feature(GenomeNode *gn, void *data, Env *env)
{
  ExtractFeatVisitor *v = (ExtractFeatVisitor*) data;
  GenomeFeatureType gf_type;
  const char *raw_sequence;
  GenomeFeature *gf;
  Range range;

  env_error_check(env);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  gf_type = genome_feature_get_type(gf);

  if (gf_type == v->type) {
    raw_sequence = bioseq_get_raw_sequence(v->bioseq);
    range = genome_node_get_range(gn);
    assert(range.start); /* 1-based coordinates */
    raw_sequence += range.start - 1;
    assert(range.end <= bioseq_get_raw_sequence_length(v->bioseq));
    str_append_cstr_nt(v->sequence, raw_sequence, range_length(range), env);
    if (genome_feature_get_strand(gf) == STRAND_REVERSE)
      v->reverse_strand = true;
  }
  return 0;
}

static int extract_feature(GenomeNode *gn, void *data, Env *env)
{
  ExtractFeatVisitor *v = (ExtractFeatVisitor*) data;
  GenomeFeatureType gf_type;
  GenomeFeature *gf;
  Range range;
  int has_err = 0;

  env_error_check(env);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  gf_type = genome_feature_get_type(gf);

  /* construct description if necessary */
  if (!str_length(v->description)) {
    str_append_cstr(v->description, genome_feature_type_get_cstr(v->type), env);
    str_append_char(v->description, '_', env);
    v->fastaseq_counter++;
    str_append_ulong(v->description, v->fastaseq_counter, env);
    if (v->join)
      str_append_cstr(v->description, " (joined)", env);
    if (v->translate)
      str_append_cstr(v->description, " (translated)", env);
  }

  if (v->join) {
    /* in this case we have to traverse the children */
    str_reset(v->sequence);
    v->reverse_strand = false;
    has_err = genome_node_traverse_direct_children(gn, v, extract_join_feature,
                                                   env);
    if (!has_err && str_length(v->sequence)) {
      if (v->reverse_strand) {
        has_err = reverse_complement(str_get(v->sequence),
                                     str_length(v->sequence), env);
      }
      if (!has_err) {
        if (v->translate) {
          str_reset(v->protein);
          translate_dna(v->protein, str_get(v->sequence),
                        str_length(v->sequence), 0, env);
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
                       range.start - 1, range_length(range), env);
    if (genome_feature_get_strand(gf) == STRAND_REVERSE) {
      has_err = reverse_complement(str_get(v->sequence),
                                   str_length(v->sequence), env);
    }
    if (!has_err) {
      if (v->translate) {
        str_reset(v->protein);
        translate_dna(v->protein, str_get(v->sequence), str_length(v->sequence),
                      0, env);
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
                                              GenomeFeature *gf, Env *env)
{
  ExtractFeatVisitor *efv;
  int has_err = 0;
  env_error_check(env);
  efv = extractfeat_visitor_cast(gv);
  if (efv->regionmapping) { /* region mapping used -> determine bioseq */
    if (!efv->sequence_file ||
        str_cmp(efv->sequence_name, genome_node_get_seqid((GenomeNode*) gf))) {
      str_delete(efv->sequence_file, env);
      efv->sequence_file = regionmapping_map(efv->regionmapping,
                              str_get(genome_node_get_seqid((GenomeNode*) gf)),
                                             env);
      if (!efv->sequence_file)
        has_err = -1;
      else {
        if (!efv->sequence_name)
          efv->sequence_name = str_new(env);
        else
          str_reset(efv->sequence_name);
        str_append_str(efv->sequence_name,
                       genome_node_get_seqid((GenomeNode*) gf), env);
        bioseq_delete(efv->bioseq, env);
        efv->bioseq = bioseq_new_str(efv->sequence_file, env);
        if (!efv->bioseq)
          has_err = -1;
      }
    }
  }
  if (!has_err) {
    has_err = genome_node_traverse_children((GenomeNode*) gf, efv,
                                            extract_feature, false, env);
  }
  return has_err;
}

static int extractfeat_visitor_sequence_region(GenomeVisitor *gv,
                                               SequenceRegion *sr, Env *env)
{
  ExtractFeatVisitor *efv;
  int has_err = 0;
  env_error_check(env);
  efv = extractfeat_visitor_cast(gv);
  if (!efv->regionmapping) { /* we have only one bioseq */
    /* check if the given sequence file contains this sequence (region) */
    if (!bioseq_contains_sequence(efv->bioseq,
                                  str_get(genome_node_get_seqid((GenomeNode*)
                                                                sr)), env)) {
      env_error_set(env,
                    "sequence \"%s\" not contained in sequence file \"%s\"",
                    str_get(genome_node_get_seqid((GenomeNode*) sr)),
                    str_get(efv->sequence_file));
      has_err = -1;
    }
  }
  return has_err;
}

const GenomeVisitorClass* extractfeat_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (ExtractFeatVisitor),
                                          extractfeat_visitor_free,
                                          NULL,
                                          extractfeat_visitor_genome_feature,
                                          extractfeat_visitor_sequence_region,
                                          NULL };
  return &gvc;
}

static GenomeVisitor* extractfeat_visitor_new(GenomeFeatureType type,
                                              bool join, bool translate,
                                              Env *env)
{
  GenomeVisitor *gv = genome_visitor_create(extractfeat_visitor_class(), env);
  ExtractFeatVisitor *efv= extractfeat_visitor_cast(gv);
  efv->description = str_new(env);
  efv->sequence = str_new(env);
  efv->protein = str_new(env);
  efv->type = type;
  efv->join = join;
  efv->translate = translate;
  efv->fastaseq_counter = 0;
  return gv;
}

GenomeVisitor* extractfeat_visitor_new_seqfile(Str *sequence_file,
                                               GenomeFeatureType type,
                                               bool join, bool translate,
                                               Env *env)
{
  GenomeVisitor *gv;
  ExtractFeatVisitor *efv;
  env_error_check(env);
  assert(sequence_file);
  gv = extractfeat_visitor_new(type, join, translate, env);
  efv = extractfeat_visitor_cast(gv);
  efv->sequence_file = str_ref(sequence_file);
  efv->bioseq = bioseq_new_str(sequence_file, env);
  if (!efv->bioseq) {
    extractfeat_visitor_free(gv, env);
    return NULL;
  }
  return gv;
}

GenomeVisitor* extractfeat_visitor_new_regionmapping(RegionMapping *rm,
                                                     GenomeFeatureType type,
                                                     bool join, bool translate,
                                                     Env *env)
{
  GenomeVisitor *gv;
  ExtractFeatVisitor *efv;
  gv = extractfeat_visitor_new(type, join, translate, env);
  efv = extractfeat_visitor_cast(gv);
  efv->regionmapping = rm;
  return gv;
}
