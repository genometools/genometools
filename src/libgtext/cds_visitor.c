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
#include "libgtcore/orf.h"
#include "libgtcore/translate.h"
#include "libgtcore/undef.h"
#include "libgtext/cds_visitor.h"
#include "libgtext/genome_visitor_rep.h"
#include "libgtext/splicedseq.h"

struct CDSVisitor {
  const GenomeVisitor parent_instance;
  Str *source;
  Splicedseq *splicedseq; /* the (spliced) sequence of the currently considered
                             gene */
  RegionMapping *region_mapping;
};

#define cds_visitor_cast(GV)\
        genome_visitor_cast(cds_visitor_class(), GV)

static void cds_visitor_free(GenomeVisitor *gv)
{
  CDSVisitor *cds_visitor = cds_visitor_cast(gv);
  assert(cds_visitor);
  str_delete(cds_visitor->source);
  splicedseq_delete(cds_visitor->splicedseq);
  region_mapping_delete(cds_visitor->region_mapping);
}

static int extract_cds_if_necessary(GenomeNode *gn, void *data, Error *err)
{
  CDSVisitor *v = (CDSVisitor*) data;
  GenomeFeature *gf;
  Range range;
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  int had_err = 0;

  error_check(err);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);

  if (genome_feature_get_type(gf) == gft_exon &&
      (genome_feature_get_strand(gf) == STRAND_FORWARD ||
       genome_feature_get_strand(gf) == STRAND_REVERSE)) {
    had_err = region_mapping_get_raw_sequence(v->region_mapping, &raw_sequence,
                                              genome_node_get_seqid(gn), err);
    if (!had_err) {
      range = genome_node_get_range(gn);
      assert(range.start && range.end); /* 1-based coordinates */
      had_err = region_mapping_get_raw_sequence_length(v->region_mapping,
                                                       &raw_sequence_length,
                                                      genome_node_get_seqid(gn),
                                                       err);
    }
    if (!had_err) {
      assert(range.end <= raw_sequence_length);
      splicedseq_add(v->splicedseq, range.start - 1, range.end - 1,
                     raw_sequence);
    }
  }
  return had_err;
}

static int extract_spliced_seq(GenomeNode *gn, CDSVisitor *visitor, Error *err)
{
  error_check(err);
  assert(gn && visitor);
  /* traverse the direct children */
  splicedseq_reset(visitor->splicedseq);
  return genome_node_traverse_direct_children(gn, visitor,
                                              extract_cds_if_necessary, err);
}

static Array* determine_ORFs_for_all_three_frames(Splicedseq *ss)
{
  Str *pr_0, *pr_1, *pr_2;
  Array *orfs;
  assert(ss);

  pr_0 = str_new();
  pr_1 = str_new();
  pr_2 = str_new();
  orfs = array_new(sizeof (Range));

  translate_dna(pr_0, splicedseq_get(ss), splicedseq_length(ss), 0);
  translate_dna(pr_1, splicedseq_get(ss), splicedseq_length(ss), 1);
  translate_dna(pr_2, splicedseq_get(ss), splicedseq_length(ss), 2);
  determine_ORFs(orfs, 0, str_get(pr_0), str_length(pr_0));
  determine_ORFs(orfs, 1, str_get(pr_1), str_length(pr_1));
  determine_ORFs(orfs, 2, str_get(pr_2), str_length(pr_2));

  str_delete(pr_2);
  str_delete(pr_1);
  str_delete(pr_0);

  return orfs;
}

static void create_CDS_features_for_ORF(Range orf, CDSVisitor *v,
                                        GenomeNode *gn)
{
  GenomeNode *cds_feature;
  unsigned long i;
  Range cds;
  Strand strand = genome_feature_get_strand((GenomeFeature*) gn);

  assert(range_length(orf) >= 3);
  /* the first CDS feature */
  cds.start = splicedseq_map(v->splicedseq, strand == STRAND_FORWARD
                             ? orf.start : orf.end) + 1;
  cds.end = splicedseq_map(v->splicedseq, strand == STRAND_FORWARD
                           ? orf.end : orf.start) + 1;
  cds_feature = genome_feature_new(gft_CDS, cds,
                                   genome_feature_get_strand((GenomeFeature*)
                                                             gn), NULL,
                                   UNDEF_ULONG);
  genome_feature_set_source(cds_feature, v->source);
  genome_node_set_seqid(cds_feature, genome_node_get_seqid(gn));
  genome_feature_set_phase(cds_feature, PHASE_ZERO);
  /* all CDS features in between */
  for (i = strand == STRAND_FORWARD ? orf.start : orf.end;
       strand == STRAND_FORWARD ? i < orf.end : i > orf.start;
       strand == STRAND_FORWARD ? i++ : i--) {
    if (splicedseq_pos_is_border(v->splicedseq, i)) {
      genome_feature_set_end((GenomeFeature*) cds_feature,
                             splicedseq_map(v->splicedseq, i) + 1);
      genome_node_is_part_of_genome_node(gn, cds_feature);
      if (strand == STRAND_FORWARD)
        orf.start = i + 1;
      else
        orf.end = i - 1;
      cds.start = splicedseq_map(v->splicedseq, strand == STRAND_FORWARD
                                 ? orf.start : orf.end) + 1;
      cds.end = splicedseq_map(v->splicedseq, strand == STRAND_FORWARD
                               ? orf.end : orf.start) + 1;
      cds_feature = genome_feature_new(gft_CDS, cds,
                                 genome_feature_get_strand((GenomeFeature*) gn),
                                       NULL, UNDEF_ULONG);
      genome_feature_set_source(cds_feature, v->source);
      genome_node_set_seqid(cds_feature, genome_node_get_seqid(gn));
      /* XXX correct this */
      genome_feature_set_phase(cds_feature, (Phase)
                               splicedseq_map(v->splicedseq, orf.start) % 3);
    }
  }
  /* set the end of the last CDS feature and store it */
  genome_feature_set_end((GenomeFeature*) cds_feature,
                         splicedseq_map(v->splicedseq,
                                        strand == STRAND_FORWARD
                                        ? orf.end : orf.start) + 1);
  genome_node_is_part_of_genome_node(gn, cds_feature);
}

static void create_CDS_features_for_longest_ORF(Array *orfs, CDSVisitor *v,
                                                GenomeNode *gn)
{
  if (array_size(orfs)) {
    /* sort ORFs according to length */
    ranges_sort_by_length_stable(orfs);

    /* create CDS features from the longest ORF */
    create_CDS_features_for_ORF(*(Range*) array_get_first(orfs), v, gn);
  }
}

static int add_cds_if_necessary(GenomeNode *gn, void *data, Error *err)
{
  CDSVisitor *v = (CDSVisitor*) data;
  GenomeFeature *gf;
  int had_err;

  error_check(err);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);

  had_err = extract_spliced_seq(gn, v, err);
  if (!had_err && splicedseq_length(v->splicedseq) > 2) {
    Array *orfs;

    if (genome_feature_get_strand(gf) == STRAND_REVERSE) {
      if (splicedseq_reverse(v->splicedseq, err))
        return -1;
    }

    orfs = determine_ORFs_for_all_three_frames(v->splicedseq);
    create_CDS_features_for_longest_ORF(orfs, v, gn);

    array_delete(orfs);
  }
  return had_err;
}

static int cds_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                      Error *err)
{
  CDSVisitor *v = cds_visitor_cast(gv);
  error_check(err);
  return genome_node_traverse_children((GenomeNode*) gf, v,
                                       add_cds_if_necessary, false, err);

}

const GenomeVisitorClass* cds_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (CDSVisitor),
                                          cds_visitor_free,
                                          NULL,
                                          cds_visitor_genome_feature,
                                          NULL };
  return &gvc;
}

GenomeVisitor* cds_visitor_new(RegionMapping *region_mapping, Str *source)
{
  GenomeVisitor *gv;
  CDSVisitor *cds_visitor;
  assert(region_mapping);
  gv = genome_visitor_create(cds_visitor_class());
  cds_visitor = cds_visitor_cast(gv);
  cds_visitor->source = str_ref(source);
  cds_visitor->splicedseq = splicedseq_new();
  cds_visitor->region_mapping = region_mapping;
  return gv;
}
