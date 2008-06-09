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
#include <stdlib.h>
#include "libgtcore/cstr.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/undef.h"
#include "libgtcore/unused.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_feature_type.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtext/genome_node_rep.h"

struct GenomeFeature
{
  const GenomeNode parent_instance;
  Str *seqid,
      *source;
  GenomeFeatureType type;
  Range range;
  double score;
  Strand strand;
  Phase phase;
  Hashtable *attributes; /* stores additional the attributes besides 'Parent'
                            and 'ID'; created on demand */
  TranscriptFeatureType transcripttype;
};

typedef struct {
  Array *exon_features,
        *cds_features;
} SaveExonAndCDSInfo;

#define genome_feature_cast(GN)\
        genome_node_cast(genome_feature_class(), GN)

static void genome_feature_free(GenomeNode *gn)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  assert(gf);
  str_delete(gf->seqid);
  str_delete(gf->source);
  hashtable_delete(gf->attributes);
}

const char* genome_feature_get_attribute(GenomeNode *gn, const char *attr_name)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  if (!gf->attributes)
    return NULL;
  return hashtable_get(gf->attributes, attr_name);
}

static Str* genome_feature_get_seqid(GenomeNode *gn)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  return gf->seqid;
}

static Range genome_feature_get_range(GenomeNode *gn)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  return gf->range;
}

static void genome_feature_set_seqid(GenomeNode *gn, Str *seqid)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  assert(gf && seqid);
  str_delete(gf->seqid);
  gf->seqid = str_ref(seqid);
}

void genome_feature_set_source(GenomeNode *gn, Str *source)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  assert(gf && source && !gf->source);
  gf->source = str_ref(source);
}

void genome_feature_set_phase(GenomeNode *gn, Phase phase)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  assert(gf && gf->phase == PHASE_UNDEFINED);
  gf->phase = phase;
}

static int genome_feature_accept(GenomeNode *gn, GenomeVisitor *gv, Error *e)
{
  GenomeFeature *gf;
  error_check(e);
  gf = genome_feature_cast(gn);
  return genome_visitor_visit_genome_feature(gv, gf, e);
}

const GenomeNodeClass* genome_feature_class()
{
  static const GenomeNodeClass gnc = { sizeof (GenomeFeature),
                                       genome_feature_free,
                                       genome_feature_get_seqid,
                                       genome_feature_get_seqid,
                                       genome_feature_get_range,
                                       NULL,
                                       genome_feature_set_seqid,
                                       genome_feature_accept };
  return &gnc;
}

GenomeNode* genome_feature_new(GenomeFeatureType type, Range range,
                               Strand strand, Str *filename,
                               unsigned long line_number)
{
  GenomeNode *gn;
  GenomeFeature *gf;
  assert(range.start <= range.end);
  gn = genome_node_create(genome_feature_class(), filename, line_number);
  gf = genome_feature_cast(gn);
  gf->seqid          = NULL;
  gf->source         = NULL;
  gf->type           = type;
  gf->score          = UNDEF_DOUBLE;
  gf->range          = range;
  gf->strand         = strand;
  gf->phase          = PHASE_UNDEFINED;
  gf->attributes     = NULL;
  gf->transcripttype = TRANSCRIPT_FEATURE_TYPE_UNDETERMINED;
  return gn;
}

GenomeNode* genome_feature_new_standard_gene(void)
{
  GenomeNode *gn, *child, *grandchild;
  Range range;
  Str *seqid;
  seqid = str_new_cstr("ctg123");

  /* gene */
  range.start = 1000; range.end = 9000;
  gn = genome_feature_new(gft_gene, range, STRAND_FORWARD, NULL, UNDEF_ULONG);
  genome_node_set_seqid(gn, seqid);

  /* TF binding site */
  range.start = 1000; range.end = 1012;
  child = genome_feature_new(gft_TF_binding_site, range, STRAND_FORWARD, NULL,
                             UNDEF_ULONG);
  genome_node_set_seqid(child, seqid);
  genome_node_is_part_of_genome_node(gn, child);

  /* first mRNA */
  range.start = 1050; range.end = 9000;
  child = genome_feature_new(gft_mRNA, range, STRAND_FORWARD, NULL,
                             UNDEF_ULONG);
  genome_node_set_seqid(child, seqid);
  genome_node_is_part_of_genome_node(gn, child);

  range.start = 1050; range.end = 1500;
  grandchild = genome_feature_new(gft_exon, range, STRAND_FORWARD, NULL,
                                  UNDEF_ULONG);
  genome_node_set_seqid(grandchild, seqid);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 3000; range.end = 3902;
  grandchild = genome_feature_new(gft_exon, range, STRAND_FORWARD, NULL,
                                  UNDEF_ULONG);
  genome_node_set_seqid(grandchild, seqid);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 5000; range.end = 5500;
  grandchild = genome_feature_new(gft_exon, range, STRAND_FORWARD, NULL,
                                  UNDEF_ULONG);
  genome_node_set_seqid(grandchild, seqid);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 7000; range.end = 9000;
  grandchild = genome_feature_new(gft_exon, range, STRAND_FORWARD, NULL,
                                  UNDEF_ULONG);
  genome_node_set_seqid(grandchild, seqid);
  genome_node_is_part_of_genome_node(child, grandchild);

  /* second mRNA */
  range.start = 1050; range.end = 9000;
  child = genome_feature_new(gft_mRNA, range, STRAND_FORWARD, NULL,
                             UNDEF_ULONG);
  genome_node_set_seqid(child, seqid);
  genome_node_is_part_of_genome_node(gn, child);

  range.start = 1050; range.end = 1500;
  grandchild = genome_feature_new(gft_exon, range, STRAND_FORWARD, NULL,
                                  UNDEF_ULONG);
  genome_node_set_seqid(grandchild, seqid);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 5000; range.end = 5500;
  grandchild = genome_feature_new(gft_exon, range, STRAND_FORWARD, NULL,
                                  UNDEF_ULONG);
  genome_node_set_seqid(grandchild, seqid);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 7000; range.end = 9000;
  grandchild = genome_feature_new(gft_exon, range, STRAND_FORWARD, NULL,
                                  UNDEF_ULONG);
  genome_node_set_seqid(grandchild, seqid);
  genome_node_is_part_of_genome_node(child, grandchild);

  /* third mRNA */
  range.start = 1300; range.end = 9000;
  child = genome_feature_new(gft_mRNA, range, STRAND_FORWARD, NULL,
                             UNDEF_ULONG);
  genome_node_set_seqid(child, seqid);
  genome_node_is_part_of_genome_node(gn, child);

  range.start = 1300; range.end = 1500;
  grandchild = genome_feature_new(gft_exon, range, STRAND_FORWARD, NULL,
                                  UNDEF_ULONG);
  genome_node_set_seqid(grandchild, seqid);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 3000; range.end = 3902;
  grandchild = genome_feature_new(gft_exon, range, STRAND_FORWARD, NULL,
                                  UNDEF_ULONG);
  genome_node_set_seqid(grandchild, seqid);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 5000; range.end = 5500;
  grandchild = genome_feature_new(gft_exon, range, STRAND_FORWARD, NULL,
                                  UNDEF_ULONG);
  genome_node_set_seqid(grandchild, seqid);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 7000; range.end = 9000;
  grandchild = genome_feature_new(gft_exon, range, STRAND_FORWARD, NULL,
                                  UNDEF_ULONG);
  genome_node_set_seqid(grandchild, seqid);
  genome_node_is_part_of_genome_node(child, grandchild);

  str_delete(seqid);
  return gn;
}

const char* genome_feature_get_source(GenomeFeature *gf)
{
  assert(gf);
  return gf->source ? str_get(gf->source) : ".";
}

GenomeFeatureType genome_feature_get_type(GenomeFeature *gf)
{
  assert(gf);
  return gf->type;
}

double genome_feature_get_score(GenomeFeature *gf)
{
  assert(gf);
  return gf->score;
}

Strand genome_feature_get_strand(GenomeFeature *gf)
{
  assert(gf);
  return gf->strand;
}

Phase genome_feature_get_phase(GenomeFeature *gf)
{
  assert(gf);
  return gf->phase;
}

static int save_exon(GenomeNode *gn, void *data, UNUSED Error *err)
{
  GenomeFeature *gf;
  Array *exon_features = (Array*) data;
  error_check(err);
  gf = (GenomeFeature*) gn;
  assert(gf && exon_features);
  if (genome_feature_get_type(gf) == gft_exon) {
    array_add(exon_features, gf);
  }
  return 0;
}

void genome_feature_get_exons(GenomeFeature *gf, Array *exon_features)
{
  int had_err;
  assert(gf && exon_features && !array_size(exon_features));
  had_err = genome_node_traverse_children((GenomeNode*) gf, exon_features,
                                          save_exon, false, NULL);
  assert(!had_err); /* cannot happen, because save_exon() is sane */
}

static int save_exons_and_cds(GenomeNode *gn, void *data, UNUSED Error *err)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  GenomeFeature *gf;
  error_check(err);
  gf = (GenomeFeature*) gn;
  assert(gf && info);
  if (genome_feature_get_type(gf) == gft_exon)
    array_add(info->exon_features, gf);
  else if (genome_feature_get_type(gf) == gft_CDS)
    array_add(info->cds_features, gf);
  return 0;
}

static void set_transcript_types(Array *features)
{
  GenomeFeature *gf;
  unsigned long i;
  assert(features);
  if (array_size(features)) {
    if (array_size(features) == 1) {
      gf = *(GenomeFeature**) array_get(features, 0);
      gf->transcripttype = TRANSCRIPT_FEATURE_TYPE_SINGLE;
    }
    else {
      gf = *(GenomeFeature**) array_get(features, 0);
      gf->transcripttype = TRANSCRIPT_FEATURE_TYPE_INITIAL;
      for (i = 1; i < array_size(features) - 1; i++) {
        gf = *(GenomeFeature**) array_get(features, i);
        gf->transcripttype = TRANSCRIPT_FEATURE_TYPE_INTERNAL;
      }
      gf = *(GenomeFeature**) array_get(features, array_size(features) - 1);
      gf->transcripttype = TRANSCRIPT_FEATURE_TYPE_TERMINAL;
    }
  }
}

static int determine_transcripttypes(GenomeNode *gn, void *data,
                                     UNUSED Error *err)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  int had_err;
  error_check(err);
  assert(gn && info);
  /* reset exon_features and cds_features */
  array_reset(info->exon_features);
  array_reset(info->cds_features);
  /* collect all direct children exons */
  had_err = genome_node_traverse_direct_children(gn, info, save_exons_and_cds,
                                                 NULL);
  assert(!had_err); /* cannot happen, because save_exon() is sane */
  /* set transcript feature type, if necessary */
  set_transcript_types(info->exon_features);
  set_transcript_types(info->cds_features);
  return 0;
}

void genome_feature_determine_transcripttypes(GenomeFeature *gf)
{
  SaveExonAndCDSInfo info;
  int had_err;
  assert(gf);
  info.exon_features = array_new(sizeof (GenomeFeature*));
  info.cds_features = array_new(sizeof (GenomeFeature*));
  had_err = genome_node_traverse_children((GenomeNode*) gf, &info,
                                          determine_transcripttypes, false,
                                          NULL);
  assert(!had_err); /* cannot happen, because determine_transcripttypes() is
                       sane */
  array_delete(info.exon_features);
  array_delete(info.cds_features);
}

TranscriptFeatureType genome_feature_get_transcriptfeaturetype(GenomeFeature
                                                               *gf)
{
  assert(gf);
  return gf->transcripttype;
}

void genome_feature_set_end(GenomeFeature *gf, unsigned long end)
{
  assert(gf && gf->range.start <= end);
  gf->range.end = end;
}

void genome_feature_set_score(GenomeFeature *gf, double score)
{
  assert(gf);
  gf->score = score;
}

void genome_feature_add_attribute(GenomeFeature *gf, const char *attr_name,
                                  const char *attr_value)
{
  assert(gf && attr_name && attr_value);
  if (!gf->attributes)
    gf->attributes = hashtable_new(HASH_STRING, ma_free_func, ma_free_func);
  hashtable_add(gf->attributes, cstr_dup(attr_name), cstr_dup(attr_value));
}

int genome_feature_foreach_attribute(GenomeFeature *gf,
                                     AttributeIterFunc iterfunc, void *data,
                                     Error *e)
{
  int had_err = 0;
  error_check(e);
  assert(gf && iterfunc);
  if (gf->attributes) {
    had_err = hashtable_foreach_ao(gf->attributes, (Hashiteratorfunc) iterfunc,
                                   data, e);
  }
  return had_err;
}

static bool genome_feature_has_gft(const GenomeFeature *gf,
                                   const GenomeFeatureType *gfts)
{
  GenomeNodeIterator *gni;
  GenomeNode *gn;
  bool has_splice_site = false;
  assert(gf && gfts && gfts[0] != undefined);
  gni = genome_node_iterator_new((GenomeNode*) gf);
  while ((gn = genome_node_iterator_next(gni))) {
    unsigned long i = 0;
    GenomeFeatureType gft = genome_feature_get_type((GenomeFeature*) gn);
    while (gfts[i] != undefined) {
      if (gft == gfts[i]) {
        has_splice_site = true;
        break;
      }
      i++;
    }
    if (has_splice_site)
      break;
  }
  genome_node_iterator_delete(gni);
  return has_splice_site;
}

bool genome_feature_has_CDS(const GenomeFeature *gf)
{
  static GenomeFeatureType gfts[] = { gft_CDS, undefined };
  return genome_feature_has_gft(gf, gfts);
}

bool genome_feature_has_splice_site(const GenomeFeature *gf)
{
  static GenomeFeatureType gfts[] = { gft_five_prime_splice_site,
                                      gft_three_prime_splice_site, undefined };
  return genome_feature_has_gft(gf, gfts);
}

double genome_feature_average_splice_site_prob(const GenomeFeature *gf)
{
  GenomeNodeIterator *gni;
  GenomeFeatureType gft;
  GenomeNode *gn;
  unsigned long num_of_splice_sites = 0;
  double averagessp = 0.0;
  assert(gf);
  gni = genome_node_iterator_new((GenomeNode*) gf);
  while ((gn = genome_node_iterator_next(gni))) {
    gft = genome_feature_get_type((GenomeFeature*) gn);
    if (gft == gft_five_prime_splice_site ||
        gft == gft_three_prime_splice_site) {
      averagessp += genome_feature_get_score((GenomeFeature*) gn);
      num_of_splice_sites++;
    }
  }
  genome_node_iterator_delete(gni);
  if (num_of_splice_sites)
    averagessp /= num_of_splice_sites;
  return averagessp;
}

bool genome_features_are_similar(GenomeFeature *gf_a, GenomeFeature *gf_b)
{
  assert(gf_a && gf_b);
  if (!str_cmp(genome_node_get_seqid((GenomeNode*) gf_a),
               genome_node_get_seqid((GenomeNode*) gf_b)) &&
      (genome_feature_get_type(gf_a) == genome_feature_get_type(gf_b)) &&
      (!range_compare(genome_feature_get_range((GenomeNode*) gf_a),
                      genome_feature_get_range((GenomeNode*) gf_b))) &&
      (genome_feature_get_strand(gf_a) == genome_feature_get_strand(gf_b)) &&
      (genome_feature_get_phase(gf_a) == genome_feature_get_phase(gf_b))) {
    return true;
  }
  return false;
}
