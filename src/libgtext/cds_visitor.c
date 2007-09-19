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
  RegionMapping *regionmapping;
};

#define cds_visitor_cast(GV)\
        genome_visitor_cast(cds_visitor_class(), GV)

static void cds_visitor_free(GenomeVisitor *gv, Env *env)
{
  CDSVisitor *cds_visitor = cds_visitor_cast(gv);
  assert(cds_visitor);
  str_delete(cds_visitor->source, env);
  splicedseq_delete(cds_visitor->splicedseq, env);
  regionmapping_delete(cds_visitor->regionmapping, env);
}

static int extract_cds_if_necessary(GenomeNode *gn, void *data, Env *env)
{
  CDSVisitor *v = (CDSVisitor*) data;
  GenomeFeature *gf;
  Range range;
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  int had_err = 0;

  env_error_check(env);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);

  if (genome_feature_get_type(gf) == gft_exon &&
      (genome_feature_get_strand(gf) == STRAND_FORWARD ||
       genome_feature_get_strand(gf) == STRAND_REVERSE)) {
    had_err = regionmapping_get_raw_sequence(v->regionmapping, &raw_sequence,
                                             genome_node_get_seqid(gn), env);
    if (!had_err) {
      range = genome_node_get_range(gn);
      assert(range.start && range.end); /* 1-based coordinates */
      had_err = regionmapping_get_raw_sequence_length(v->regionmapping,
                                                      &raw_sequence_length,
                                                      genome_node_get_seqid(gn),
                                                      env);
    }
    if (!had_err) {
      assert(range.end <= raw_sequence_length);
      splicedseq_add(v->splicedseq, range.start - 1, range.end - 1,
                     raw_sequence, env);
    }
  }
  return had_err;
}

static int add_cds_if_necessary(GenomeNode *gn, void *data, Env *env)
{
  CDSVisitor *v = (CDSVisitor*) data;
  GenomeNode *cds_feature;
  Str *pr_0, *pr_1, *pr_2;
  GenomeFeature *gf;
  unsigned long i;
  Array *orfs;
  Range orf, cds;
  Strand strand;
  int had_err;

  env_error_check(env);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);

  /* traverse the direct children */
  splicedseq_reset(v->splicedseq);
  had_err = genome_node_traverse_direct_children(gn, v,
                                                 extract_cds_if_necessary, env);
  if (!had_err && splicedseq_length(v->splicedseq) > 2) {
    strand = genome_feature_get_strand(gf);
    if (strand == STRAND_REVERSE) {
      if (splicedseq_reverse(v->splicedseq, env))
        return -1;
    }
    /* determine ORFs for all three frames */
    pr_0 = str_new(env);
    pr_1 = str_new(env);
    pr_2 = str_new(env);
    /* printf("pr_0=%s\n", str_get(pr_0)); */
    orfs = array_new(sizeof (Range), env);
    translate_dna(pr_0, splicedseq_get(v->splicedseq),
                  splicedseq_length(v->splicedseq), 0, env);
    translate_dna(pr_1, splicedseq_get(v->splicedseq),
                  splicedseq_length(v->splicedseq), 1, env);
    translate_dna(pr_2, splicedseq_get(v->splicedseq),
                  splicedseq_length(v->splicedseq), 2, env);
    determine_ORFs(orfs, 0, str_get(pr_0), str_length(pr_0), env);
    determine_ORFs(orfs, 1, str_get(pr_1), str_length(pr_1), env);
    determine_ORFs(orfs, 2, str_get(pr_2), str_length(pr_2), env);

    if (array_size(orfs)) {
      /* sort ORFs according to length */
      ranges_sort_by_length_stable(orfs, env);

      /* create CDS features from the longest ORF */
      orf = *(Range*) array_get(orfs, 0);
      assert(range_length(orf) >= 3);
      /* the first CDS feature */
      /*printf("%lu, %lu\n", orf.start, orf.end); */
      cds.start = splicedseq_map(v->splicedseq, strand == STRAND_FORWARD
                                 ? orf.start : orf.end) + 1;
      cds.end = splicedseq_map(v->splicedseq, strand == STRAND_FORWARD
                               ? orf.end : orf.start) + 1;
      /*printf("%lu, %lu\n", cds.start, cds.end);*/
      cds_feature = genome_feature_new(gft_CDS, cds,
                                       genome_feature_get_strand(gf), NULL,
                                       UNDEF_ULONG, env);
      genome_node_set_source(cds_feature, v->source);
      genome_node_set_seqid(cds_feature, genome_node_get_seqid(gn), env);
      genome_node_set_phase(cds_feature, PHASE_ZERO);
      /* all CDS features in between */
      for (i = strand == STRAND_FORWARD ? orf.start + 1 : orf.end - 1;
           strand == STRAND_FORWARD ? i < orf.end : i > orf.start;
           strand == STRAND_FORWARD ? i++ : i--) {
        if (splicedseq_pos_is_border(v->splicedseq, i)) {
          /*printf("i=%lu\n", i);*/
          genome_feature_set_end((GenomeFeature*) cds_feature,
                                 splicedseq_map(v->splicedseq, i) + 1);
          genome_node_is_part_of_genome_node(gn, cds_feature, env);
          if (strand == STRAND_FORWARD)
            orf.start = i + 1;
          else
            orf.end = i - 1;
          cds.start = splicedseq_map(v->splicedseq, strand == STRAND_FORWARD
                                     ? orf.start : orf.end) + 1;
          cds.end = splicedseq_map(v->splicedseq, strand == STRAND_FORWARD
                                   ? orf.end : orf.start) + 1;
          cds_feature = genome_feature_new(gft_CDS, cds,
                                           genome_feature_get_strand(gf),
                                           NULL, UNDEF_ULONG, env);
          genome_node_set_source(cds_feature, v->source);
          genome_node_set_seqid(cds_feature, genome_node_get_seqid(gn), env);
          /* XXX correct this */
          genome_node_set_phase(cds_feature, (Phase)
                                splicedseq_map(v->splicedseq, orf.start) % 3);
        }
      }
      /* set the end of the last CDS feature and store it */
      genome_feature_set_end((GenomeFeature*) cds_feature,
                             splicedseq_map(v->splicedseq,
                                            strand == STRAND_FORWARD
                                            ? orf.end : orf.start) + 1);
      genome_node_is_part_of_genome_node(gn, cds_feature, env);
    }

    /* free */
    array_delete(orfs, env);
    str_delete(pr_2, env);
    str_delete(pr_1, env);
    str_delete(pr_0, env);
  }
  return had_err;
}

static int cds_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                      Env *env)
{
  CDSVisitor *v = cds_visitor_cast(gv);
  env_error_check(env);
  return genome_node_traverse_children((GenomeNode*) gf, v,
                                       add_cds_if_necessary, false, env);

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

GenomeVisitor* cds_visitor_new(RegionMapping *regionmapping, Str *source,
                               Env *env)
{
  GenomeVisitor *gv;
  CDSVisitor *cds_visitor;
  env_error_check(env);
  assert(regionmapping);
  gv = genome_visitor_create(cds_visitor_class(), env);
  cds_visitor = cds_visitor_cast(gv);
  cds_visitor->source = str_ref(source);
  cds_visitor->splicedseq = splicedseq_new(env);
  cds_visitor->regionmapping = regionmapping;
  return gv;
}
