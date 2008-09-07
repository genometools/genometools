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

#define GT_STRAND_OFFSET                   5
#define GT_STRAND_MASK                     0x7
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

struct GT_GenomeFeature
{
  const GT_GenomeNode parent_instance;
  GT_Str *seqid,
      *source;
  GT_FeatureTypeFactory *ftf;
  GT_GenomeFeatureType *type;
  GT_Range range;
  float score;
  TagValueMap attributes; /* stores the attributes; created on demand */
  GT_GenomeFeature *representative;
};

typedef struct {
  GT_Array *exon_features,
        *cds_features;
} SaveExonAndCDSInfo;

#define gt_genome_feature_cast(GN)\
        gt_genome_node_cast(gt_genome_feature_class(), GN)

static void gt_genome_feature_free(GT_GenomeNode *gn)
{
  GT_GenomeFeature *gf = gt_genome_feature_cast(gn);
  assert(gf);
  str_delete(gf->seqid);
  str_delete(gf->source);
  gt_feature_type_factory_delete(gf->ftf);
  tag_value_map_delete(gf->attributes);
}

const char* gt_genome_feature_get_attribute(GT_GenomeNode *gn, const char *attr_name)
{
  GT_GenomeFeature *gf = gt_genome_feature_cast(gn);
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

GT_StrArray* gt_genome_feature_get_attribute_list(GT_GenomeFeature *gf)
{
  GT_StrArray *list = gt_strarray_new();
  if (gf->attributes)
    tag_value_map_foreach(gf->attributes, store_attribute, list);
  return list;
}

static GT_Str* gt_genome_feature_get_seqid(GT_GenomeNode *gn)
{
  GT_GenomeFeature *gf = gt_genome_feature_cast(gn);
  return gf->seqid;
}

static GT_Range gt_genome_feature_get_range(GT_GenomeNode *gn)
{
  GT_GenomeFeature *gf = gt_genome_feature_cast(gn);
  return gf->range;
}

static void gt_genome_feature_set_range(GT_GenomeNode *gn, GT_Range range)
{
  GT_GenomeFeature *gf = gt_genome_feature_cast(gn);
  gf->range = range;
}

static void gt_genome_feature_change_seqid(GT_GenomeNode *gn, GT_Str *seqid)
{
  GT_GenomeFeature *gf = gt_genome_feature_cast(gn);
  assert(gf && seqid);
  str_delete(gf->seqid);
  gf->seqid = str_ref(seqid);
}

void gt_genome_feature_set_source(GT_GenomeNode *gn, GT_Str *source)
{
  GT_GenomeFeature *gf = gt_genome_feature_cast(gn);
  assert(gf && source && !gf->source);
  gf->source = str_ref(source);
}

void gt_genome_feature_set_phase(GT_GenomeNode *gn, Phase phase)
{
  assert(gn);
  gn->bit_field &= ~(PHASE_MASK << PHASE_OFFSET);
  gn->bit_field |= phase << PHASE_OFFSET;
}

static int gt_genome_feature_accept(GT_GenomeNode *gn, GenomeVisitor *gv, GT_Error *err)
{
  GT_GenomeFeature *gf;
  gt_error_check(err);
  gf = gt_genome_feature_cast(gn);
  return genome_visitor_visit_genome_feature(gv, gf, err);
}

const GT_GenomeNodeClass* gt_genome_feature_class()
{
  static const GT_GenomeNodeClass gnc = { sizeof (GT_GenomeFeature),
                                       gt_genome_feature_free,
                                       gt_genome_feature_get_seqid,
                                       gt_genome_feature_get_seqid,
                                       gt_genome_feature_get_range,
                                       gt_genome_feature_set_range,
                                       gt_genome_feature_change_seqid,
                                       gt_genome_feature_accept };
  return &gnc;
}

static void set_transcriptfeaturetype(GT_GenomeNode *gn, TranscriptFeatureType tft)
{
  assert(gn);
  gn->bit_field &= ~(TRANSCRIPT_FEATURE_TYPE_MASK <<
                     TRANSCRIPT_FEATURE_TYPE_OFFSET);
  gn->bit_field |= tft << TRANSCRIPT_FEATURE_TYPE_OFFSET;
}

GT_GenomeNode* gt_genome_feature_new(GT_Str *seqid, GT_GenomeFeatureType *type, GT_Range range,
                               GT_Strand strand)
{
  GT_GenomeNode *gn;
  GT_GenomeFeature *gf;
  assert(seqid && type);
  assert(range.start <= range.end);
  gn = gt_genome_node_create(gt_genome_feature_class());
  gf = gt_genome_feature_cast(gn);
  gf->seqid     = str_ref(seqid);
  gf->source    = NULL;
  gf->ftf       = gt_feature_type_factory_ref(gt_genome_feature_type_get_ftf(type));
  gf->type      = type;
  gf->score     = UNDEF_FLOAT;
  gf->range     = range;
  gn->bit_field |= strand << GT_STRAND_OFFSET;
  gt_genome_feature_set_phase(gn, PHASE_UNDEFINED);
  set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_UNDETERMINED);
  gf->attributes     = NULL;
  gf->representative = NULL;
  return gn;
}

GT_GenomeNode* gt_genome_feature_new_pseudo(GT_GenomeFeature *gf)
{
  GT_GenomeFeature *pf;
  GT_GenomeNode *pn;
  assert(gf);
  pn = gt_genome_feature_new(gt_genome_feature_get_seqid((GT_GenomeNode*) gf),
                          gt_genome_feature_get_type(gf),
                          gt_genome_feature_get_range((GT_GenomeNode*) gf),
                          gt_genome_feature_get_strand(gf));
  pf = gt_genome_feature_cast(pn);
  pf->type = NULL; /* pseudo features do not have a type */
  gt_genome_feature_set_source(pn, gf->source);
  pn->bit_field |= 1 << PSEUDO_FEATURE_OFFSET;
  return pn;
}

GT_GenomeNode* gt_genome_feature_new_standard_gene(GT_FeatureTypeFactory *ftf)
{
  GT_GenomeNode *gn, *child, *grandchild;
  GT_GenomeFeatureType *type;
  GT_Range range;
  GT_Str *seqid;
  seqid = str_new_cstr("ctg123");

  /* gene */
  range.start = 1000; range.end = 9000;
  type = gt_feature_type_factory_create_gft(ftf, gft_gene);
  gn = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);

  /* TF binding site */
  range.start = 1000; range.end = 1012;
  type = gt_feature_type_factory_create_gft(ftf, gft_TF_binding_site);
  child = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(gn, child);

  /* first mRNA */
  range.start = 1050; range.end = 9000;
  type = gt_feature_type_factory_create_gft(ftf, gft_mRNA);
  child = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(gn, child);

  range.start = 1050; range.end = 1500;
  type = gt_feature_type_factory_create_gft(ftf, gft_exon);
  grandchild = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 3000; range.end = 3902;
  grandchild = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 5000; range.end = 5500;
  grandchild = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 7000; range.end = 9000;
  grandchild = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(child, grandchild);

  /* second mRNA */
  range.start = 1050; range.end = 9000;
  type = gt_feature_type_factory_create_gft(ftf, gft_mRNA);
  child = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(gn, child);

  range.start = 1050; range.end = 1500;
  type = gt_feature_type_factory_create_gft(ftf, gft_exon);
  grandchild = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 5000; range.end = 5500;
  grandchild = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 7000; range.end = 9000;
  grandchild = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(child, grandchild);

  /* third mRNA */
  range.start = 1300; range.end = 9000;
  type = gt_feature_type_factory_create_gft(ftf, gft_mRNA);
  child = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(gn, child);

  range.start = 1300; range.end = 1500;
  type = gt_feature_type_factory_create_gft(ftf, gft_exon);
  grandchild = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 3000; range.end = 3902;
  grandchild = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 5000; range.end = 5500;
  grandchild = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(child, grandchild);

  range.start = 7000; range.end = 9000;
  grandchild = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  gt_genome_node_is_part_of_genome_node(child, grandchild);

  str_delete(seqid);
  return gn;
}

const char* gt_genome_feature_get_source(GT_GenomeFeature *gf)
{
  assert(gf);
  return gf->source ? str_get(gf->source) : ".";
}

GT_GenomeFeatureType* gt_genome_feature_get_type(GT_GenomeFeature *gf)
{
  assert(gf);
  return gf->type;
}

GT_GenomeFeatureType* gt_genome_feature_create_gft(GT_GenomeFeature *gf,
                                             const char *type)
{
  assert(gf && type);
  return gt_genome_feature_type_create_gft(gf->type, type);
}

bool gt_genome_feature_has_type(GT_GenomeFeature *gf, const char *type)
{
  assert(gf && type);
  return gt_genome_feature_type_is(gf->type, type);
}

bool gt_genome_feature_score_is_defined(const GT_GenomeFeature *gf)
{
  GT_GenomeNode *gn = (GT_GenomeNode*) gf;
  assert(gn);
  if ((gn->bit_field >> SCORE_IS_DEFINED_OFFSET) & SCORE_IS_DEFINED_MASK)
    return true;
  return false;
}

bool gt_genome_feature_is_multi(const GT_GenomeFeature *gf)
{
  GT_GenomeNode *gn = (GT_GenomeNode*) gf;
  assert(gn);
  if ((gn->bit_field >> MULTI_FEATURE_OFFSET) & MULTI_FEATURE_MASK) {
    assert(!gt_genome_feature_is_pseudo(gf));
    return true;
  }
  return false;
}

bool gt_genome_feature_is_pseudo(const GT_GenomeFeature *gf)
{
  GT_GenomeNode *gn = (GT_GenomeNode*) gf;
  assert(gn);
  if ((gn->bit_field >> PSEUDO_FEATURE_OFFSET) & PSEUDO_FEATURE_MASK) {
    assert(!gt_genome_feature_is_multi(gf));
    return true;
  }
  return false;
}

static void gt_genome_feature_set_multi(const GT_GenomeFeature *gf)
{
  GT_GenomeNode *gn;
  assert(gf && !gt_genome_feature_is_multi(gf));
  gn = (GT_GenomeNode*) gf;
  gn->bit_field |= 1 << MULTI_FEATURE_OFFSET;
}

void gt_genome_feature_make_multi_representative(const GT_GenomeFeature *gf)
{
  assert(gf && !gt_genome_feature_is_multi(gf));
  gt_genome_feature_set_multi(gf);
}

void gt_genome_feature_set_multi_representative(GT_GenomeFeature *gf,
                                             GT_GenomeFeature *rep)
{
  assert(gf && !gt_genome_feature_is_multi(gf));
  assert(rep && gt_genome_feature_is_multi(rep));
  gt_genome_feature_set_multi(gf);
  gf->representative = rep;
}

GT_GenomeFeature* gt_genome_feature_get_multi_representative(GT_GenomeFeature *gf)
{
  assert(gf && gt_genome_feature_is_multi(gf) && !gt_genome_feature_is_pseudo(gf));
  if (gf->representative) {
    assert(gt_genome_feature_is_multi(gf->representative));
    assert(gt_genome_feature_get_multi_representative(gf->representative) ==
           gf->representative);
    return gf->representative;
  }
  return gf; /* is itself the representative */
}

float gt_genome_feature_get_score(GT_GenomeFeature *gf)
{
  assert(gf);
  assert(gt_genome_feature_score_is_defined(gf));
  return gf->score;
}

GT_Strand gt_genome_feature_get_strand(GT_GenomeFeature *gf)
{
  GT_GenomeNode *gn = (GT_GenomeNode*) gf;
  assert(gf);
  return (gn->bit_field >> GT_STRAND_OFFSET) & GT_STRAND_MASK;
}

Phase gt_genome_feature_get_phase(GT_GenomeFeature *gf)
{
  GT_GenomeNode *gn = (GT_GenomeNode*) gf;
  assert(gn);
  return (gn->bit_field >> PHASE_OFFSET) & PHASE_MASK;
}

static int save_exon(GT_GenomeNode *gn, void *data, UNUSED GT_Error *err)
{
  GT_GenomeFeature *gf;
  GT_Array *exon_features = (GT_Array*) data;
  gt_error_check(err);
  gf = (GT_GenomeFeature*) gn;
  assert(gf && exon_features);
  if (gt_genome_feature_has_type(gf, gft_exon)) {
    gt_array_add(exon_features, gf);
  }
  return 0;
}

void gt_genome_feature_get_exons(GT_GenomeFeature *gf, GT_Array *exon_features)
{
  int had_err;
  assert(gf && exon_features && !gt_array_size(exon_features));
  had_err = gt_genome_node_traverse_children((GT_GenomeNode*) gf, exon_features,
                                          save_exon, false, NULL);
  assert(!had_err); /* cannot happen, because save_exon() is sane */
}

static int save_exons_and_cds(GT_GenomeNode *gn, void *data, UNUSED GT_Error *err)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  GT_GenomeFeature *gf;
  gt_error_check(err);
  gf = (GT_GenomeFeature*) gn;
  assert(gf && info);
  if (gt_genome_feature_has_type(gf, gft_exon))
    gt_array_add(info->exon_features, gf);
  else if (gt_genome_feature_has_type(gf, gft_CDS))
    gt_array_add(info->cds_features, gf);
  return 0;
}

static void set_transcript_types(GT_Array *features)
{
  GT_GenomeNode *gn;
  unsigned long i;
  assert(features);
  if (gt_array_size(features)) {
    if (gt_array_size(features) == 1) {
      gn = *(GT_GenomeNode**) gt_array_get(features, 0);
      set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_SINGLE);
    }
    else {
      gn = *(GT_GenomeNode**) gt_array_get(features, 0);
      set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_INITIAL);
      for (i = 1; i < gt_array_size(features) - 1; i++) {
        gn = *(GT_GenomeNode**) gt_array_get(features, i);
        set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_INTERNAL);
      }
      gn = *(GT_GenomeNode**) gt_array_get(features, gt_array_size(features) - 1);
      set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_TERMINAL);
    }
  }
}

static int determine_transcripttypes(GT_GenomeNode *gn, void *data,
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
  had_err = gt_genome_node_traverse_direct_children(gn, info, save_exons_and_cds,
                                                 NULL);
  assert(!had_err); /* cannot happen, because save_exon() is sane */
  /* set transcript feature type, if necessary */
  set_transcript_types(info->exon_features);
  set_transcript_types(info->cds_features);
  return 0;
}

void gt_genome_feature_determine_transcripttypes(GT_GenomeFeature *gf)
{
  SaveExonAndCDSInfo info;
  int had_err;
  assert(gf);
  info.exon_features = gt_array_new(sizeof (GT_GenomeFeature*));
  info.cds_features = gt_array_new(sizeof (GT_GenomeFeature*));
  had_err = gt_genome_node_traverse_children((GT_GenomeNode*) gf, &info,
                                          determine_transcripttypes, false,
                                          NULL);
  assert(!had_err); /* cannot happen, because determine_transcripttypes() is
                       sane */
  gt_array_delete(info.exon_features);
  gt_array_delete(info.cds_features);
}

TranscriptFeatureType gt_genome_feature_get_transcriptfeaturetype(GT_GenomeFeature
                                                               *gf)
{
  GT_GenomeNode *gn = (GT_GenomeNode*) gf;
  assert(gn);
  return (gn->bit_field >> TRANSCRIPT_FEATURE_TYPE_OFFSET) &
         TRANSCRIPT_FEATURE_TYPE_MASK;
}

void gt_genome_feature_set_end(GT_GenomeFeature *gf, unsigned long end)
{
  assert(gf && gf->range.start <= end);
  gf->range.end = end;
}

void gt_genome_feature_set_score(GT_GenomeFeature *gf, float score)
{
  GT_GenomeNode *gn = (GT_GenomeNode*) gf;
  assert(gf);
  gn->bit_field |= 1 << SCORE_IS_DEFINED_OFFSET;
  gf->score = score;
}

void gt_genome_feature_unset_score(GT_GenomeFeature *gf)
{
  GT_GenomeNode *gn = (GT_GenomeNode*) gf;
  assert(gf);
  gn->bit_field &= ~(1 << SCORE_IS_DEFINED_OFFSET);
  gf->score = UNDEF_FLOAT;
}

void gt_genome_feature_add_attribute(GT_GenomeFeature *gf, const char *attr_name,
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

void gt_genome_feature_foreach_attribute(GT_GenomeFeature *gf,
                                      AttributeIterFunc iterfunc, void *data)
{
  assert(gf && iterfunc);
  if (gf->attributes) {
    tag_value_map_foreach(gf->attributes, (TagValueMapIteratorFunc) iterfunc,
                          data);
  }
}

static bool gt_genome_feature_has_gft(const GT_GenomeFeature *gf, const char **gfts)
{
  GT_GenomeNodeIterator *gni;
  GT_GenomeNode *gn;
  bool has_gft = false;
  assert(gf && gfts && gfts[0] != NULL);
  gni = gt_genome_node_iterator_new((GT_GenomeNode*) gf);
  while ((gn = gt_genome_node_iterator_next(gni))) {
    unsigned long i = 0;
    while (gfts[i] != NULL) {
      if (gt_genome_feature_has_type((GT_GenomeFeature*) gn, gfts[i])) {
        has_gft = true;
        break;
      }
      i++;
    }
    if (has_gft)
      break;
  }
  gt_genome_node_iterator_delete(gni);
  return has_gft;
}

bool gt_genome_feature_has_CDS(const GT_GenomeFeature *gf)
{
  static const char *gfts[] = { gft_CDS, NULL };
  return gt_genome_feature_has_gft(gf, gfts);
}

bool gt_genome_feature_has_splice_site(const GT_GenomeFeature *gf)
{
  static const char *gfts[] = { gft_five_prime_splice_site,
                                gft_three_prime_splice_site, NULL };
  return gt_genome_feature_has_gft(gf, gfts);
}

double gt_genome_feature_average_splice_site_prob(const GT_GenomeFeature *gf)
{
  GT_GenomeNodeIterator *gni;
  GT_GenomeNode *gn;
  unsigned long num_of_splice_sites = 0;
  double averagessp = 0.0;
  assert(gf);
  gni = gt_genome_node_iterator_new((GT_GenomeNode*) gf);
  while ((gn = gt_genome_node_iterator_next(gni))) {
    if (gt_genome_feature_has_type((GT_GenomeFeature*) gn,
                                gft_five_prime_splice_site) ||
        gt_genome_feature_has_type((GT_GenomeFeature*) gn,
                                gft_three_prime_splice_site)) {
      averagessp += gt_genome_feature_get_score((GT_GenomeFeature*) gn);
      num_of_splice_sites++;
    }
  }
  gt_genome_node_iterator_delete(gni);
  if (num_of_splice_sites)
    averagessp /= num_of_splice_sites;
  return averagessp;
}

bool genome_features_are_similar(GT_GenomeFeature *gf_a, GT_GenomeFeature *gf_b)
{
  assert(gf_a && gf_b);
  if (!str_cmp(gt_genome_node_get_seqid((GT_GenomeNode*) gf_a),
               gt_genome_node_get_seqid((GT_GenomeNode*) gf_b)) &&
      (gt_genome_feature_get_type(gf_a) == gt_genome_feature_get_type(gf_b)) &&
      (!gt_range_compare(gt_genome_feature_get_range((GT_GenomeNode*) gf_a),
                      gt_genome_feature_get_range((GT_GenomeNode*) gf_b))) &&
      (gt_genome_feature_get_strand(gf_a) == gt_genome_feature_get_strand(gf_b)) &&
      (gt_genome_feature_get_phase(gf_a) == gt_genome_feature_get_phase(gf_b))) {
    return true;
  }
  return false;
}

int gt_genome_feature_unit_test(GT_Error *err)
{
  GT_FeatureTypeFactory *feature_type_factory;
  GT_GenomeFeatureType *type;
  GT_GenomeNode *gf;
  GT_Range range;
  GT_Str *seqid;
  int had_err = 0;

  gt_error_check(err);

  seqid = str_new_cstr("seqid");
  feature_type_factory = gt_feature_type_factory_builtin_new();
  type = gt_feature_type_factory_create_gft(feature_type_factory, "gene");
  range.start = 1;
  range.end = 1000;
  gf = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);

  ensure(had_err, !gt_genome_feature_score_is_defined((GT_GenomeFeature*) gf));

  gt_genome_node_delete(gf);
  gt_feature_type_factory_delete(feature_type_factory);
  str_delete(seqid);

  return had_err;
}
