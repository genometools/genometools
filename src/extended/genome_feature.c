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
#include <string.h>
#include "core/cstr.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/undef.h"
#include "core/unused.h"
#include "extended/feature_type_factory.h"
#include "extended/feature_type_factory_builtin.h"
#include "extended/genome_feature.h"
#include "extended/genome_feature_type.h"
#include "extended/genome_feature_type_imp.h"
#include "extended/genome_node_iterator.h"
#include "extended/genome_node_rep.h"
#include "extended/tag_value_map.h"

#define STRAND_OFFSET                   5
#define STRAND_MASK                     0x7
#define PHASE_OFFSET                    8
#define PHASE_MASK                      0x3
#define TRANSCRIPT_FEATURE_TYPE_OFFSET  10
#define TRANSCRIPT_FEATURE_TYPE_MASK    0x7
#define SCORE_IS_DEFINED_OFFSET         13
#define SCORE_IS_DEFINED_MASK           0x1
#define MULTI_FEATURE_OFFSET            14
#define MULTI_FEATURE_MASK              0x1
#define PSEUDO_FEATURE_OFFSET           15
#define PSEUDO_FEATURE_MASK             0x1

struct GenomeFeature
{
  const GenomeNode parent_instance;
  Str *seqid,
      *source;
  FeatureTypeFactory *ftf;
  GenomeFeatureType *type;
  Range range;
  float score;
  TagValueMap attributes; /* stores the attributes; created on demand */
  GenomeFeature *representative;
};

typedef struct {
  GT_Array *exon_features,
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
  feature_type_factory_delete(gf->ftf);
  tag_value_map_delete(gf->attributes);
}

const char* genome_feature_get_attribute(GenomeNode *gn, const char *attr_name)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  if (!gf->attributes)
    return NULL;
  return tag_value_map_get(gf->attributes, attr_name);
}

static void store_attribute(const char *attr_name,
                            UNUSED const char *attr_value, void *data)
{
  GT_StrArray *list = data;
  assert(attr_name && attr_value && data);
  gt_strarray_add_cstr(list, attr_name);
}

GT_StrArray* genome_feature_get_attribute_list(GenomeFeature *gf)
{
  GT_StrArray *list = gt_strarray_new();
  if (gf->attributes)
    tag_value_map_foreach(gf->attributes, store_attribute, list);
  return list;
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

static void genome_feature_set_range(GenomeNode *gn, Range range)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  gf->range = range;
}

static void genome_feature_change_seqid(GenomeNode *gn, Str *seqid)
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
  assert(gn);
  gn->bit_field &= ~(PHASE_MASK << PHASE_OFFSET);
  gn->bit_field |= phase << PHASE_OFFSET;
}

static int genome_feature_accept(GenomeNode *gn, GenomeVisitor *gv, GT_Error *err)
{
  GenomeFeature *gf;
  gt_error_check(err);
  gf = genome_feature_cast(gn);
  return genome_visitor_visit_genome_feature(gv, gf, err);
}

const GenomeNodeClass* genome_feature_class()
{
  static const GenomeNodeClass gnc = { sizeof (GenomeFeature),
                                       genome_feature_free,
                                       genome_feature_get_seqid,
                                       genome_feature_get_seqid,
                                       genome_feature_get_range,
                                       genome_feature_set_range,
                                       genome_feature_change_seqid,
                                       genome_feature_accept };
  return &gnc;
}

static void set_transcriptfeaturetype(GenomeNode *gn, TranscriptFeatureType tft)
{
  assert(gn);
  gn->bit_field &= ~(TRANSCRIPT_FEATURE_TYPE_MASK <<
                     TRANSCRIPT_FEATURE_TYPE_OFFSET);
  gn->bit_field |= tft << TRANSCRIPT_FEATURE_TYPE_OFFSET;
}

GenomeNode* genome_feature_new(Str *seqid, GenomeFeatureType *type, Range range,
                               Strand strand)
{
  GenomeNode *gn;
  GenomeFeature *gf;
  assert(seqid && type);
  assert(range.start <= range.end);
  gn = genome_node_create(genome_feature_class());
  gf = genome_feature_cast(gn);
  gf->seqid     = str_ref(seqid);
  gf->source    = NULL;
  gf->ftf       = feature_type_factory_ref(genome_feature_type_get_ftf(type));
  gf->type      = type;
  gf->score     = UNDEF_FLOAT;
  gf->range     = range;
  gn->bit_field |= strand << STRAND_OFFSET;
  genome_feature_set_phase(gn, PHASE_UNDEFINED);
  set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_UNDETERMINED);
  gf->attributes     = NULL;
  gf->representative = NULL;
  return gn;
}

GenomeNode* genome_feature_new_pseudo(GenomeFeature *gf)
{
  GenomeFeature *pf;
  GenomeNode *pn;
  assert(gf);
  pn = genome_feature_new(genome_feature_get_seqid((GenomeNode*) gf),
                          genome_feature_get_type(gf),
                          genome_feature_get_range((GenomeNode*) gf),
                          genome_feature_get_strand(gf));
  pf = genome_feature_cast(pn);
  pf->type = NULL; /* pseudo features do not have a type */
  genome_feature_set_source(pn, gf->source);
  pn->bit_field |= 1 << PSEUDO_FEATURE_OFFSET;
  return pn;
}

GenomeNode* genome_feature_new_standard_gene(FeatureTypeFactory *ftf)
{
  GenomeNode *gn, *child, *grandchild;
  GenomeFeatureType *type;
  Range range;
  Str *seqid;
  seqid = str_new_cstr("ctg123");

  /* gene */
  range.start = 1000; range.end = 9000;
  type = feature_type_factory_create_gft(ftf, gft_gene);
  gn = genome_feature_new(seqid, type, range, STRAND_FORWARD);

  /* TF binding site */
  range.start = 1000; range.end = 1012;
  type = feature_type_factory_create_gft(ftf, gft_TF_binding_site);
  child = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(gn, child);

  /* first mRNA */
  range.start = 1050; range.end = 9000;
  type = feature_type_factory_create_gft(ftf, gft_mRNA);
  child = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(gn, child);

  range.start = 1050; range.end = 1500;
  type = feature_type_factory_create_gft(ftf, gft_exon);
  grandchild = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 3000; range.end = 3902;
  grandchild = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 5000; range.end = 5500;
  grandchild = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 7000; range.end = 9000;
  grandchild = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(child, grandchild);

  /* second mRNA */
  range.start = 1050; range.end = 9000;
  type = feature_type_factory_create_gft(ftf, gft_mRNA);
  child = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(gn, child);

  range.start = 1050; range.end = 1500;
  type = feature_type_factory_create_gft(ftf, gft_exon);
  grandchild = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 5000; range.end = 5500;
  grandchild = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 7000; range.end = 9000;
  grandchild = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(child, grandchild);

  /* third mRNA */
  range.start = 1300; range.end = 9000;
  type = feature_type_factory_create_gft(ftf, gft_mRNA);
  child = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(gn, child);

  range.start = 1300; range.end = 1500;
  type = feature_type_factory_create_gft(ftf, gft_exon);
  grandchild = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 3000; range.end = 3902;
  grandchild = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 5000; range.end = 5500;
  grandchild = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 7000; range.end = 9000;
  grandchild = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  genome_node_is_part_of_genome_node(child, grandchild);

  str_delete(seqid);
  return gn;
}

const char* genome_feature_get_source(GenomeFeature *gf)
{
  assert(gf);
  return gf->source ? str_get(gf->source) : ".";
}

GenomeFeatureType* genome_feature_get_type(GenomeFeature *gf)
{
  assert(gf);
  return gf->type;
}

GenomeFeatureType* genome_feature_create_gft(GenomeFeature *gf,
                                             const char *type)
{
  assert(gf && type);
  return genome_feature_type_create_gft(gf->type, type);
}

bool genome_feature_has_type(GenomeFeature *gf, const char *type)
{
  assert(gf && type);
  return genome_feature_type_is(gf->type, type);
}

bool genome_feature_score_is_defined(const GenomeFeature *gf)
{
  GenomeNode *gn = (GenomeNode*) gf;
  assert(gn);
  if ((gn->bit_field >> SCORE_IS_DEFINED_OFFSET) & SCORE_IS_DEFINED_MASK)
    return true;
  return false;
}

bool genome_feature_is_multi(const GenomeFeature *gf)
{
  GenomeNode *gn = (GenomeNode*) gf;
  assert(gn);
  if ((gn->bit_field >> MULTI_FEATURE_OFFSET) & MULTI_FEATURE_MASK) {
    assert(!genome_feature_is_pseudo(gf));
    return true;
  }
  return false;
}

bool genome_feature_is_pseudo(const GenomeFeature *gf)
{
  GenomeNode *gn = (GenomeNode*) gf;
  assert(gn);
  if ((gn->bit_field >> PSEUDO_FEATURE_OFFSET) & PSEUDO_FEATURE_MASK) {
    assert(!genome_feature_is_multi(gf));
    return true;
  }
  return false;
}

static void genome_feature_set_multi(const GenomeFeature *gf)
{
  GenomeNode *gn;
  assert(gf && !genome_feature_is_multi(gf));
  gn = (GenomeNode*) gf;
  gn->bit_field |= 1 << MULTI_FEATURE_OFFSET;
}

void genome_feature_make_multi_representative(const GenomeFeature *gf)
{
  assert(gf && !genome_feature_is_multi(gf));
  genome_feature_set_multi(gf);
}

void genome_feature_set_multi_representative(GenomeFeature *gf,
                                             GenomeFeature *rep)
{
  assert(gf && !genome_feature_is_multi(gf));
  assert(rep && genome_feature_is_multi(rep));
  genome_feature_set_multi(gf);
  gf->representative = rep;
}

GenomeFeature* genome_feature_get_multi_representative(GenomeFeature *gf)
{
  assert(gf && genome_feature_is_multi(gf) && !genome_feature_is_pseudo(gf));
  if (gf->representative) {
    assert(genome_feature_is_multi(gf->representative));
    assert(genome_feature_get_multi_representative(gf->representative) ==
           gf->representative);
    return gf->representative;
  }
  return gf; /* is itself the representative */
}

float genome_feature_get_score(GenomeFeature *gf)
{
  assert(gf);
  assert(genome_feature_score_is_defined(gf));
  return gf->score;
}

Strand genome_feature_get_strand(GenomeFeature *gf)
{
  GenomeNode *gn = (GenomeNode*) gf;
  assert(gf);
  return (gn->bit_field >> STRAND_OFFSET) & STRAND_MASK;
}

Phase genome_feature_get_phase(GenomeFeature *gf)
{
  GenomeNode *gn = (GenomeNode*) gf;
  assert(gn);
  return (gn->bit_field >> PHASE_OFFSET) & PHASE_MASK;
}

static int save_exon(GenomeNode *gn, void *data, UNUSED GT_Error *err)
{
  GenomeFeature *gf;
  GT_Array *exon_features = (GT_Array*) data;
  gt_error_check(err);
  gf = (GenomeFeature*) gn;
  assert(gf && exon_features);
  if (genome_feature_has_type(gf, gft_exon)) {
    gt_array_add(exon_features, gf);
  }
  return 0;
}

void genome_feature_get_exons(GenomeFeature *gf, GT_Array *exon_features)
{
  int had_err;
  assert(gf && exon_features && !gt_array_size(exon_features));
  had_err = genome_node_traverse_children((GenomeNode*) gf, exon_features,
                                          save_exon, false, NULL);
  assert(!had_err); /* cannot happen, because save_exon() is sane */
}

static int save_exons_and_cds(GenomeNode *gn, void *data, UNUSED GT_Error *err)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  GenomeFeature *gf;
  gt_error_check(err);
  gf = (GenomeFeature*) gn;
  assert(gf && info);
  if (genome_feature_has_type(gf, gft_exon))
    gt_array_add(info->exon_features, gf);
  else if (genome_feature_has_type(gf, gft_CDS))
    gt_array_add(info->cds_features, gf);
  return 0;
}

static void set_transcript_types(GT_Array *features)
{
  GenomeNode *gn;
  unsigned long i;
  assert(features);
  if (gt_array_size(features)) {
    if (gt_array_size(features) == 1) {
      gn = *(GenomeNode**) gt_array_get(features, 0);
      set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_SINGLE);
    }
    else {
      gn = *(GenomeNode**) gt_array_get(features, 0);
      set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_INITIAL);
      for (i = 1; i < gt_array_size(features) - 1; i++) {
        gn = *(GenomeNode**) gt_array_get(features, i);
        set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_INTERNAL);
      }
      gn = *(GenomeNode**) gt_array_get(features, gt_array_size(features) - 1);
      set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_TERMINAL);
    }
  }
}

static int determine_transcripttypes(GenomeNode *gn, void *data,
                                     UNUSED GT_Error *err)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  int had_err;
  gt_error_check(err);
  assert(gn && info);
  /* reset exon_features and cds_features */
  gt_array_reset(info->exon_features);
  gt_array_reset(info->cds_features);
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
  info.exon_features = gt_array_new(sizeof (GenomeFeature*));
  info.cds_features = gt_array_new(sizeof (GenomeFeature*));
  had_err = genome_node_traverse_children((GenomeNode*) gf, &info,
                                          determine_transcripttypes, false,
                                          NULL);
  assert(!had_err); /* cannot happen, because determine_transcripttypes() is
                       sane */
  gt_array_delete(info.exon_features);
  gt_array_delete(info.cds_features);
}

TranscriptFeatureType genome_feature_get_transcriptfeaturetype(GenomeFeature
                                                               *gf)
{
  GenomeNode *gn = (GenomeNode*) gf;
  assert(gn);
  return (gn->bit_field >> TRANSCRIPT_FEATURE_TYPE_OFFSET) &
         TRANSCRIPT_FEATURE_TYPE_MASK;
}

void genome_feature_set_end(GenomeFeature *gf, unsigned long end)
{
  assert(gf && gf->range.start <= end);
  gf->range.end = end;
}

void genome_feature_set_score(GenomeFeature *gf, float score)
{
  GenomeNode *gn = (GenomeNode*) gf;
  assert(gf);
  gn->bit_field |= 1 << SCORE_IS_DEFINED_OFFSET;
  gf->score = score;
}

void genome_feature_unset_score(GenomeFeature *gf)
{
  GenomeNode *gn = (GenomeNode*) gf;
  assert(gf);
  gn->bit_field &= ~(1 << SCORE_IS_DEFINED_OFFSET);
  gf->score = UNDEF_FLOAT;
}

void genome_feature_add_attribute(GenomeFeature *gf, const char *attr_name,
                                  const char *attr_value)
{
  assert(gf && attr_name && attr_value);
  assert(strlen(attr_name)); /* attribute name cannot be empty */
  assert(strlen(attr_value)); /* attribute value cannot be empty */
  if (!gf->attributes)
    gf->attributes = tag_value_map_new(attr_name, attr_value);
  else
    tag_value_map_add(&gf->attributes, attr_name, attr_value);
}

void genome_feature_foreach_attribute(GenomeFeature *gf,
                                      AttributeIterFunc iterfunc, void *data)
{
  assert(gf && iterfunc);
  if (gf->attributes) {
    tag_value_map_foreach(gf->attributes, (TagValueMapIteratorFunc) iterfunc,
                          data);
  }
}

static bool genome_feature_has_gft(const GenomeFeature *gf, const char **gfts)
{
  GenomeNodeIterator *gni;
  GenomeNode *gn;
  bool has_gft = false;
  assert(gf && gfts && gfts[0] != NULL);
  gni = genome_node_iterator_new((GenomeNode*) gf);
  while ((gn = genome_node_iterator_next(gni))) {
    unsigned long i = 0;
    while (gfts[i] != NULL) {
      if (genome_feature_has_type((GenomeFeature*) gn, gfts[i])) {
        has_gft = true;
        break;
      }
      i++;
    }
    if (has_gft)
      break;
  }
  genome_node_iterator_delete(gni);
  return has_gft;
}

bool genome_feature_has_CDS(const GenomeFeature *gf)
{
  static const char *gfts[] = { gft_CDS, NULL };
  return genome_feature_has_gft(gf, gfts);
}

bool genome_feature_has_splice_site(const GenomeFeature *gf)
{
  static const char *gfts[] = { gft_five_prime_splice_site,
                                gft_three_prime_splice_site, NULL };
  return genome_feature_has_gft(gf, gfts);
}

double genome_feature_average_splice_site_prob(const GenomeFeature *gf)
{
  GenomeNodeIterator *gni;
  GenomeNode *gn;
  unsigned long num_of_splice_sites = 0;
  double averagessp = 0.0;
  assert(gf);
  gni = genome_node_iterator_new((GenomeNode*) gf);
  while ((gn = genome_node_iterator_next(gni))) {
    if (genome_feature_has_type((GenomeFeature*) gn,
                                gft_five_prime_splice_site) ||
        genome_feature_has_type((GenomeFeature*) gn,
                                gft_three_prime_splice_site)) {
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

int genome_feature_unit_test(GT_Error *err)
{
  FeatureTypeFactory *feature_type_factory;
  GenomeFeatureType *type;
  GenomeNode *gf;
  Range range;
  Str *seqid;
  int had_err = 0;

  gt_error_check(err);

  seqid = str_new_cstr("seqid");
  feature_type_factory = feature_type_factory_builtin_new();
  type = feature_type_factory_create_gft(feature_type_factory, "gene");
  range.start = 1;
  range.end = 1000;
  gf = genome_feature_new(seqid, type, range, STRAND_FORWARD);

  ensure(had_err, !genome_feature_score_is_defined((GenomeFeature*) gf));

  genome_node_delete(gf);
  feature_type_factory_delete(feature_type_factory);
  str_delete(seqid);

  return had_err;
}
