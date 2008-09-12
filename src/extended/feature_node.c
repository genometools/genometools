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
#include "core/strcmp.h"
#include "core/symbol.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
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

struct GtFeatureNode
{
  const GtGenomeNode parent_instance;
  GtStr *seqid,
         *source;
  const char *type;
  GtRange range;
  float score;
  TagValueMap attributes; /* stores the attributes; created on demand */
  GtFeatureNode *representative;
};

typedef struct {
  GtArray *exon_features,
        *cds_features;
} SaveExonAndCDSInfo;

static void gt_feature_node_free(GtGenomeNode *gn)
{
  GtFeatureNode *gf = gt_feature_node_cast(gn);
  assert(gf);
  gt_str_delete(gf->seqid);
  gt_str_delete(gf->source);
  tag_value_map_delete(gf->attributes);
}

const char* gt_feature_node_get_attribute(GtGenomeNode *gn,
                                            const char *attr_name)
{
  GtFeatureNode *gf = gt_feature_node_cast(gn);
  if (!gf->attributes)
    return NULL;
  return tag_value_map_get(gf->attributes, attr_name);
}

static void store_attribute(const char *attr_name,
                            GT_UNUSED const char *attr_value, void *data)
{
  GtStrArray *list = data;
  assert(attr_name && attr_value && data);
  gt_strarray_add_cstr(list, attr_name);
}

GtStrArray* gt_feature_node_get_attribute_list(GtFeatureNode *gf)
{
  GtStrArray *list = gt_strarray_new();
  if (gf->attributes)
    tag_value_map_foreach(gf->attributes, store_attribute, list);
  return list;
}

static GtStr* gt_feature_node_get_seqid(GtGenomeNode *gn)
{
  GtFeatureNode *gf = gt_feature_node_cast(gn);
  return gf->seqid;
}

static GtRange gt_feature_node_get_range(GtGenomeNode *gn)
{
  GtFeatureNode *gf = gt_feature_node_cast(gn);
  return gf->range;
}

static void gt_feature_node_set_range(GtGenomeNode *gn, GtRange range)
{
  GtFeatureNode *gf = gt_feature_node_cast(gn);
  gf->range = range;
}

static void gt_feature_node_change_seqid(GtGenomeNode *gn, GtStr *seqid)
{
  GtFeatureNode *gf = gt_feature_node_cast(gn);
  assert(gf && seqid);
  gt_str_delete(gf->seqid);
  gf->seqid = gt_str_ref(seqid);
}

void gt_feature_node_set_source(GtGenomeNode *gn, GtStr *source)
{
  GtFeatureNode *gf = gt_feature_node_cast(gn);
  assert(gf && source && !gf->source);
  gf->source = gt_str_ref(source);
}

void gt_feature_node_set_phase(GtGenomeNode *gn, Phase phase)
{
  assert(gn);
  gn->bit_field &= ~(PHASE_MASK << PHASE_OFFSET);
  gn->bit_field |= phase << PHASE_OFFSET;
}

static int gt_feature_node_accept(GtGenomeNode *gn, GenomeVisitor *gv,
                                  GtError *err)
{
  GtFeatureNode *gf;
  gt_error_check(err);
  gf = gt_feature_node_cast(gn);
  return genome_visitor_visit_feature_node(gv, gf, err);
}

const GtGenomeNodeClass* gt_feature_node_class()
{
  static const GtGenomeNodeClass gnc = { sizeof (GtFeatureNode),
                                       gt_feature_node_free,
                                       gt_feature_node_get_seqid,
                                       gt_feature_node_get_seqid,
                                       gt_feature_node_get_range,
                                       gt_feature_node_set_range,
                                       gt_feature_node_change_seqid,
                                       gt_feature_node_accept };
  return &gnc;
}

static void set_transcriptfeaturetype(GtGenomeNode *gn,
                                      TranscriptFeatureType tft)
{
  assert(gn);
  gn->bit_field &= ~(TRANSCRIPT_FEATURE_TYPE_MASK <<
                     TRANSCRIPT_FEATURE_TYPE_OFFSET);
  gn->bit_field |= tft << TRANSCRIPT_FEATURE_TYPE_OFFSET;
}

GtGenomeNode* gt_feature_node_new(GtStr *seqid, const char *type,
                                     unsigned long start, unsigned long end,
                                     GtStrand strand)
{
  GtGenomeNode *gn;
  GtFeatureNode *gf;
  assert(seqid && type);
  assert(start <= end);
  gn = gt_genome_node_create(gt_feature_node_class());
  gf = gt_feature_node_cast(gn);
  gf->seqid       = gt_str_ref(seqid);
  gf->source      = NULL;
  gf->type        = gt_symbol(type);
  gf->score       = UNDEF_FLOAT;
  gf->range.start = start;
  gf->range.end   = end;
  gn->bit_field |= strand << GT_STRAND_OFFSET;
  gt_feature_node_set_phase(gn, GT_PHASE_UNDEFINED);
  set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_UNDETERMINED);
  gf->attributes     = NULL;
  gf->representative = NULL;
  return gn;
}

GtGenomeNode* gt_feature_node_new_pseudo(GtFeatureNode *gf)
{
  GtFeatureNode *pf;
  GtGenomeNode *pn;
  GtRange range;
  assert(gf);
  range = gt_feature_node_get_range((GtGenomeNode*) gf),
  pn = gt_feature_node_new(gt_feature_node_get_seqid((GtGenomeNode*) gf),
                            gt_feature_node_get_type(gf), range.start,
                            range.end, gt_feature_node_get_strand(gf));
  pf = gt_feature_node_cast(pn);
  pf->type = NULL; /* pseudo features do not have a type */
  gt_feature_node_set_source(pn, gf->source);
  pn->bit_field |= 1 << PSEUDO_FEATURE_OFFSET;
  return pn;
}

GtGenomeNode* gt_feature_node_new_standard_gene(void)
{
  GtGenomeNode *gn, *child, *grandchild;
  GtStr *seqid;
  seqid = gt_str_new_cstr("ctg123");

  /* gene */
  gn = gt_feature_node_new(seqid, gft_gene, 1000, 9000, GT_STRAND_FORWARD);

  /* TF binding site */
  child = gt_feature_node_new(seqid, gft_TF_binding_site, 1000, 1012,
                                GT_STRAND_FORWARD);
  gt_genome_node_add_child(gn, child);

  /* first mRNA */
  child = gt_feature_node_new(seqid, gft_mRNA, 1050, 9000, GT_STRAND_FORWARD);
  gt_genome_node_add_child(gn, child);

  grandchild = gt_feature_node_new(seqid, gft_exon, 1050, 1500,
                                     GT_STRAND_FORWARD);
  gt_genome_node_add_child(child, grandchild);

  grandchild = gt_feature_node_new(seqid, gft_exon, 3000, 3902,
                                     GT_STRAND_FORWARD);
  gt_genome_node_add_child(child, grandchild);

  grandchild = gt_feature_node_new(seqid, gft_exon, 5000, 5500,
                                     GT_STRAND_FORWARD);
  gt_genome_node_add_child(child, grandchild);

  grandchild = gt_feature_node_new(seqid, gft_exon, 7000, 9000,
                                     GT_STRAND_FORWARD);
  gt_genome_node_add_child(child, grandchild);

  /* second mRNA */
  child = gt_feature_node_new(seqid, gft_mRNA, 1050, 9000, GT_STRAND_FORWARD);
  gt_genome_node_add_child(gn, child);

  grandchild = gt_feature_node_new(seqid, gft_exon, 1050, 1500,
                                     GT_STRAND_FORWARD);
  gt_genome_node_add_child(child, grandchild);

  grandchild = gt_feature_node_new(seqid, gft_exon, 5000, 5500,
                                     GT_STRAND_FORWARD);
  gt_genome_node_add_child(child, grandchild);

  grandchild = gt_feature_node_new(seqid, gft_exon, 7000, 9000,
                                     GT_STRAND_FORWARD);
  gt_genome_node_add_child(child, grandchild);

  /* third mRNA */
  child = gt_feature_node_new(seqid, gft_mRNA, 1300, 9000, GT_STRAND_FORWARD);
  gt_genome_node_add_child(gn, child);

  grandchild = gt_feature_node_new(seqid, gft_exon, 1300, 1500,
                                     GT_STRAND_FORWARD);
  gt_genome_node_add_child(child, grandchild);

  grandchild = gt_feature_node_new(seqid, gft_exon, 3000, 3902,
                                     GT_STRAND_FORWARD);
  gt_genome_node_add_child(child, grandchild);

  grandchild = gt_feature_node_new(seqid, gft_exon, 5000, 5500,
                                     GT_STRAND_FORWARD);
  gt_genome_node_add_child(child, grandchild);

  grandchild = gt_feature_node_new(seqid, gft_exon, 7000, 9000,
                                     GT_STRAND_FORWARD);
  gt_genome_node_add_child(child, grandchild);

  gt_str_delete(seqid);
  return gn;
}

const char* gt_feature_node_get_source(GtFeatureNode *gf)
{
  assert(gf);
  return gf->source ? gt_str_get(gf->source) : ".";
}

const char* gt_feature_node_get_type(GtFeatureNode *gf)
{
  assert(gf);
  return gf->type;
}

bool gt_feature_node_has_type(GtFeatureNode *gf, const char *type)
{
  assert(gf && type);
  return gt_strcmp(gf->type, type) ? false : true;
}

bool gt_feature_node_score_is_defined(const GtFeatureNode *gf)
{
  GtGenomeNode *gn = (GtGenomeNode*) gf;
  assert(gn);
  if ((gn->bit_field >> SCORE_IS_DEFINED_OFFSET) & SCORE_IS_DEFINED_MASK)
    return true;
  return false;
}

bool gt_feature_node_is_multi(const GtFeatureNode *gf)
{
  GtGenomeNode *gn = (GtGenomeNode*) gf;
  assert(gn);
  if ((gn->bit_field >> MULTI_FEATURE_OFFSET) & MULTI_FEATURE_MASK) {
    assert(!gt_feature_node_is_pseudo(gf));
    return true;
  }
  return false;
}

bool gt_feature_node_is_pseudo(const GtFeatureNode *gf)
{
  GtGenomeNode *gn = (GtGenomeNode*) gf;
  assert(gn);
  if ((gn->bit_field >> PSEUDO_FEATURE_OFFSET) & PSEUDO_FEATURE_MASK) {
    assert(!gt_feature_node_is_multi(gf));
    return true;
  }
  return false;
}

static void gt_feature_node_set_multi(const GtFeatureNode *gf)
{
  GtGenomeNode *gn;
  assert(gf && !gt_feature_node_is_multi(gf));
  gn = (GtGenomeNode*) gf;
  gn->bit_field |= 1 << MULTI_FEATURE_OFFSET;
}

void gt_feature_node_make_multi_representative(const GtFeatureNode *gf)
{
  assert(gf && !gt_feature_node_is_multi(gf));
  gt_feature_node_set_multi(gf);
}

void gt_feature_node_set_multi_representative(GtFeatureNode *gf,
                                             GtFeatureNode *rep)
{
  assert(gf && !gt_feature_node_is_multi(gf));
  assert(rep && gt_feature_node_is_multi(rep));
  gt_feature_node_set_multi(gf);
  gf->representative = rep;
}

GtFeatureNode* gt_feature_node_get_multi_representative(GtFeatureNode
                                                             *gf)
{
  assert(gf && gt_feature_node_is_multi(gf) &&
         !gt_feature_node_is_pseudo(gf));
  if (gf->representative) {
    assert(gt_feature_node_is_multi(gf->representative));
    assert(gt_feature_node_get_multi_representative(gf->representative) ==
           gf->representative);
    return gf->representative;
  }
  return gf; /* is itself the representative */
}

float gt_feature_node_get_score(GtFeatureNode *gf)
{
  assert(gf);
  assert(gt_feature_node_score_is_defined(gf));
  return gf->score;
}

GtStrand gt_feature_node_get_strand(GtFeatureNode *gf)
{
  GtGenomeNode *gn = (GtGenomeNode*) gf;
  assert(gf);
  return (gn->bit_field >> GT_STRAND_OFFSET) & GT_STRAND_MASK;
}

Phase gt_feature_node_get_phase(GtFeatureNode *gf)
{
  GtGenomeNode *gn = (GtGenomeNode*) gf;
  assert(gn);
  return (gn->bit_field >> PHASE_OFFSET) & PHASE_MASK;
}

static int save_exon(GtGenomeNode *gn, void *data, GT_UNUSED GtError *err)
{
  GtFeatureNode *gf;
  GtArray *exon_features = (GtArray*) data;
  gt_error_check(err);
  gf = (GtFeatureNode*) gn;
  assert(gf && exon_features);
  if (gt_feature_node_has_type(gf, gft_exon)) {
    gt_array_add(exon_features, gf);
  }
  return 0;
}

void gt_feature_node_get_exons(GtFeatureNode *gf, GtArray *exon_features)
{
  int had_err;
  assert(gf && exon_features && !gt_array_size(exon_features));
  had_err = gt_genome_node_traverse_children((GtGenomeNode*) gf, exon_features,
                                          save_exon, false, NULL);
  assert(!had_err); /* cannot happen, because save_exon() is sane */
}

static int save_exons_and_cds(GtGenomeNode *gn, void *data,
                              GT_UNUSED GtError *err)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  GtFeatureNode *gf;
  gt_error_check(err);
  gf = (GtFeatureNode*) gn;
  assert(gf && info);
  if (gt_feature_node_has_type(gf, gft_exon))
    gt_array_add(info->exon_features, gf);
  else if (gt_feature_node_has_type(gf, gft_CDS))
    gt_array_add(info->cds_features, gf);
  return 0;
}

static void set_transcript_types(GtArray *features)
{
  GtGenomeNode *gn;
  unsigned long i;
  assert(features);
  if (gt_array_size(features)) {
    if (gt_array_size(features) == 1) {
      gn = *(GtGenomeNode**) gt_array_get(features, 0);
      set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_SINGLE);
    }
    else {
      gn = *(GtGenomeNode**) gt_array_get(features, 0);
      set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_INITIAL);
      for (i = 1; i < gt_array_size(features) - 1; i++) {
        gn = *(GtGenomeNode**) gt_array_get(features, i);
        set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_INTERNAL);
      }
      gn = *(GtGenomeNode**)
           gt_array_get(features, gt_array_size(features) - 1);
      set_transcriptfeaturetype(gn, TRANSCRIPT_FEATURE_TYPE_TERMINAL);
    }
  }
}

static int determine_transcripttypes(GtGenomeNode *gn, void *data,
                                     GT_UNUSED GtError *err)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  int had_err;
  gt_error_check(err);
  assert(gn && info);
  /* reset exon_features and cds_features */
  gt_array_reset(info->exon_features);
  gt_array_reset(info->cds_features);
  /* collect all direct children exons */
  had_err = gt_genome_node_traverse_direct_children(gn, info,
                                                    save_exons_and_cds, NULL);
  assert(!had_err); /* cannot happen, because save_exons_and_cds() is sane */
  /* set transcript feature type, if necessary */
  set_transcript_types(info->exon_features);
  set_transcript_types(info->cds_features);
  return 0;
}

void gt_feature_node_determine_transcripttypes(GtFeatureNode *gf)
{
  SaveExonAndCDSInfo info;
  int had_err;
  assert(gf);
  info.exon_features = gt_array_new(sizeof (GtFeatureNode*));
  info.cds_features = gt_array_new(sizeof (GtFeatureNode*));
  had_err = gt_genome_node_traverse_children((GtGenomeNode*) gf, &info,
                                          determine_transcripttypes, false,
                                          NULL);
  assert(!had_err); /* cannot happen, because determine_transcripttypes() is
                       sane */
  gt_array_delete(info.exon_features);
  gt_array_delete(info.cds_features);
}

TranscriptFeatureType gt_feature_node_get_transcriptfeaturetype(
                                                           GtFeatureNode *gf)
{
  GtGenomeNode *gn = (GtGenomeNode*) gf;
  assert(gn);
  return (gn->bit_field >> TRANSCRIPT_FEATURE_TYPE_OFFSET) &
         TRANSCRIPT_FEATURE_TYPE_MASK;
}

void gt_feature_node_set_end(GtFeatureNode *gf, unsigned long end)
{
  assert(gf && gf->range.start <= end);
  gf->range.end = end;
}

void gt_feature_node_set_score(GtFeatureNode *gf, float score)
{
  GtGenomeNode *gn = (GtGenomeNode*) gf;
  assert(gf);
  gn->bit_field |= 1 << SCORE_IS_DEFINED_OFFSET;
  gf->score = score;
}

void gt_feature_node_unset_score(GtFeatureNode *gf)
{
  GtGenomeNode *gn = (GtGenomeNode*) gf;
  assert(gf);
  gn->bit_field &= ~(1 << SCORE_IS_DEFINED_OFFSET);
  gf->score = UNDEF_FLOAT;
}

void gt_feature_node_add_attribute(GtFeatureNode *gf,
                                     const char *attr_name,
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

void gt_feature_node_foreach_attribute(GtFeatureNode *gf,
                                      AttributeIterFunc iterfunc, void *data)
{
  assert(gf && iterfunc);
  if (gf->attributes) {
    tag_value_map_foreach(gf->attributes, (TagValueMapIteratorFunc) iterfunc,
                          data);
  }
}

static bool genome_feature_has_gft(const GtFeatureNode *gf,
                                   const char **gfts)
{
  GtGenomeNodeIterator *gni;
  GtGenomeNode *gn;
  bool has_gft = false;
  assert(gf && gfts && gfts[0] != NULL);
  gni = gt_genome_node_iterator_new((GtGenomeNode*) gf);
  while ((gn = gt_genome_node_iterator_next(gni))) {
    unsigned long i = 0;
    while (gfts[i] != NULL) {
      if (gt_feature_node_has_type((GtFeatureNode*) gn, gfts[i])) {
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

bool gt_feature_node_has_CDS(const GtFeatureNode *gf)
{
  static const char *gfts[] = { gft_CDS, NULL };
  return genome_feature_has_gft(gf, gfts);
}

bool gt_feature_node_has_splice_site(const GtFeatureNode *gf)
{
  static const char *gfts[] = { gft_five_prime_splice_site,
                                gft_three_prime_splice_site, NULL };
  return genome_feature_has_gft(gf, gfts);
}

double gt_feature_node_average_splice_site_prob(const GtFeatureNode *gf)
{
  GtGenomeNodeIterator *gni;
  GtGenomeNode *gn;
  unsigned long num_of_splice_sites = 0;
  double averagessp = 0.0;
  assert(gf);
  gni = gt_genome_node_iterator_new((GtGenomeNode*) gf);
  while ((gn = gt_genome_node_iterator_next(gni))) {
    if (gt_feature_node_has_type((GtFeatureNode*) gn,
                                gft_five_prime_splice_site) ||
        gt_feature_node_has_type((GtFeatureNode*) gn,
                                gft_three_prime_splice_site)) {
      averagessp += gt_feature_node_get_score((GtFeatureNode*) gn);
      num_of_splice_sites++;
    }
  }
  gt_genome_node_iterator_delete(gni);
  if (num_of_splice_sites)
    averagessp /= num_of_splice_sites;
  return averagessp;
}

bool gt_genome_features_are_similar(GtFeatureNode *gf_a,
                                    GtFeatureNode *gf_b)
{
  assert(gf_a && gf_b);
  if (!gt_str_cmp(gt_genome_node_get_seqid((GtGenomeNode*) gf_a),
                  gt_genome_node_get_seqid((GtGenomeNode*) gf_b)) &&
      (gt_feature_node_get_type(gf_a) == gt_feature_node_get_type(gf_b)) &&
      (!gt_range_compare(gt_feature_node_get_range((GtGenomeNode*) gf_a),
                         gt_feature_node_get_range((GtGenomeNode*) gf_b))) &&
      (gt_feature_node_get_strand(gf_a) ==
       gt_feature_node_get_strand(gf_b)) &&
      (gt_feature_node_get_phase(gf_a) ==
       gt_feature_node_get_phase(gf_b))) {
    return true;
  }
  return false;
}

int gt_feature_node_unit_test(GtError *err)
{
  GtGenomeNode *gf;
  GtStr *seqid;
  int had_err = 0;

  gt_error_check(err);

  seqid = gt_str_new_cstr("seqid");
  gf = gt_feature_node_new(seqid, gft_gene, 1, 1000, GT_STRAND_FORWARD);

  ensure(had_err, !gt_feature_node_score_is_defined((GtFeatureNode*) gf));

  gt_genome_node_delete(gf);
  gt_str_delete(seqid);

  return had_err;
}
