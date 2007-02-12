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

struct Extractfeat_visitor {
  const Genome_visitor parent_instance;
  Str *sequence_file,
      *description, /* the description of the currently extracted feature */
      *sequence,    /* the sequence of the currently extracted feature */
      *protein;
  Genome_feature_type type;
  bool join,
       translate,
       reverse_strand;
  Bioseq *bioseq;
  unsigned long fastaseq_counter;
};

#define extractfeat_visitor_cast(GV)\
        genome_visitor_cast(extractfeat_visitor_class(), GV)

static void extractfeat_visitor_free(Genome_visitor *gv)
{
  Extractfeat_visitor *extractfeat_visitor = extractfeat_visitor_cast(gv);
  assert(extractfeat_visitor);
  str_free(extractfeat_visitor->sequence_file);
  str_free(extractfeat_visitor->description);
  str_free(extractfeat_visitor->sequence);
  str_free(extractfeat_visitor->protein);
  bioseq_free(extractfeat_visitor->bioseq);
}

static void extract_join_feature(Genome_node *gn, void *data)
{
  Extractfeat_visitor *v = (Extractfeat_visitor*) data;
  Genome_feature_type gf_type;
  const char *raw_sequence;
  Genome_feature *gf;
  Range range;

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
}

static void extract_feature(Genome_node *gn, void *data)
{
  Extractfeat_visitor *v = (Extractfeat_visitor*) data;
  Genome_feature_type gf_type;
  Genome_feature *gf;
  Range range;

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
    genome_node_traverse_direct_children(gn, v, extract_join_feature);
    if (str_length(v->sequence)) {
      if (v->reverse_strand)
        reverse_complement(str_get(v->sequence), str_length(v->sequence));
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
      str_reset(v->description);
    }
  }
  else if (gf_type == v->type) {
    /* otherwise we only have to look this feature */
    range = genome_node_get_range(gn);
    assert(range.start); /* 1-based coordinates */
    assert(range.end <= bioseq_get_raw_sequence_length(v->bioseq));
    str_reset(v->sequence);
    str_append_cstr_nt(v->sequence, bioseq_get_raw_sequence(v->bioseq) +
                       range.start - 1, range_length(range));
    if (genome_feature_get_strand(gf) == STRAND_REVERSE)
      reverse_complement(str_get(v->sequence), str_length(v->sequence));
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
    str_reset(v->description);
  }
}

static void extractfeat_visitor_genome_feature(Genome_visitor *gv,
                                               Genome_feature *gf,
                                               /*@unused@*/ Log *l)
{
  Extractfeat_visitor *v = extractfeat_visitor_cast(gv);
  genome_node_traverse_children((Genome_node*) gf, v, extract_feature, 0);
}

static void extractfeat_visitor_sequence_region(Genome_visitor *gv,
                                                Sequence_region *sr,
                                                /*@unused@*/ Log *l)
{
  Extractfeat_visitor *extractfeat_visitor = extractfeat_visitor_cast(gv);
  /* check if the given sequence file contains this sequence (region) */
  if (!bioseq_contains_sequence(extractfeat_visitor->bioseq,
                                str_get(genome_node_get_seqid((Genome_node*)
                                                              sr)))) {
    error("sequence \"%s\" not contained in sequence file \"%s\"",
          str_get(genome_node_get_seqid((Genome_node*) sr)),
          str_get(extractfeat_visitor->sequence_file));
  }
}

const Genome_visitor_class* extractfeat_visitor_class()
{
  static const Genome_visitor_class gvc = { sizeof(Extractfeat_visitor),
                                            extractfeat_visitor_free,
                                            NULL,
                                            extractfeat_visitor_genome_feature,
                                            extractfeat_visitor_sequence_region,
                                            NULL };
  return &gvc;
}

Genome_visitor* extractfeat_visitor_new(Str *sequence_file,
                                        Genome_feature_type type,
                                        bool join, bool translate)
{
  Genome_visitor *gv = genome_visitor_create(extractfeat_visitor_class());
  Extractfeat_visitor *extractfeat_visitor = extractfeat_visitor_cast(gv);
  extractfeat_visitor->sequence_file = str_ref(sequence_file);
  extractfeat_visitor->description = str_new();
  extractfeat_visitor->sequence = str_new();
  extractfeat_visitor->protein = str_new();
  extractfeat_visitor->type = type;
  extractfeat_visitor->join = join;
  extractfeat_visitor->translate = translate;
  extractfeat_visitor->bioseq = bioseq_new_str(sequence_file);
  extractfeat_visitor->fastaseq_counter = 0;
  return gv;
}
