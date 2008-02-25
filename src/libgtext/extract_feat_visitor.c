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
#include "libgtext/extract_feat_sequence.h"
#include "libgtext/extract_feat_visitor.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtext/genome_visitor_rep.h"

struct ExtractFeatVisitor {
  const GenomeVisitor parent_instance;
  GenomeFeatureType type;
  bool join,
       translate;
  unsigned long fastaseq_counter;
  RegionMapping *region_mapping;
};

#define extract_feat_visitor_cast(GV)\
        genome_visitor_cast(extract_feat_visitor_class(), GV)

static void extract_feat_visitor_free(GenomeVisitor *gv)
{
  ExtractFeatVisitor *extract_feat_visitor = extract_feat_visitor_cast(gv);
  assert(extract_feat_visitor);
  region_mapping_delete(extract_feat_visitor->region_mapping);
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

static void show_entry(Str *description, Str *sequence, bool translate)
{
  if (translate) {
    Str *protein = str_new();
    translate_dna(protein, str_get(sequence), str_length(sequence), 0);
    fasta_show_entry(str_get(description), str_get(protein),
                     str_length(protein), 0);
    str_delete(protein);
  }
  else {
    fasta_show_entry(str_get(description), str_get(sequence),
                     str_length(sequence), 0);
  }
}

static int extract_feat_visitor_genome_feature(GenomeVisitor *gv,
                                               GenomeFeature *gf, Error *err)
{
  ExtractFeatVisitor *efv;
  GenomeNodeIterator *gni;
  GenomeNode *gn;
  Str *description,
      *sequence;
  int had_err = 0;
  error_check(err);
  efv = extract_feat_visitor_cast(gv);
  assert(efv->region_mapping);
  gni = genome_node_iterator_new((GenomeNode*) gf);
  description = str_new();
  sequence = str_new();
  while (!had_err && (gn = genome_node_iterator_next(gni))) {
    if (extract_feat_sequence(sequence, gn, efv->type, efv->join,
                              efv->region_mapping, err)) {
      had_err = -1;
    }

    if (!had_err && str_length(sequence)) {
      efv->fastaseq_counter++;
      construct_description(description, efv->type, efv->fastaseq_counter,
                            efv->join, efv->translate);
      show_entry(description, sequence, efv->translate);
      str_reset(description);
      str_reset(sequence);
    }
  }
  str_delete(sequence);
  str_delete(description);
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
  efv->type = type;
  efv->join = join;
  efv->translate = translate;
  efv->fastaseq_counter = 0;
  efv->region_mapping = rm;
  return gv;
}
