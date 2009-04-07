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

#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/cstr.h"
#include "core/ensure.h"
#include "core/hashtable.h"
#include "core/ma.h"
#include "core/queue.h"
#include "core/strcmp.h"
#include "core/symbol.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_rep.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/genome_node_rep.h"
#include "extended/tag_value_map.h"

#define PARENT_STATUS_OFFSET            1
#define PARENT_STATUS_MASK              0x3
#define TREE_STATUS_OFFSET              3
#define TREE_STATUS_MASK                0x3
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

typedef enum {
  NO_PARENT,
  ONE_PARENT,
  MULTIPLE_PARENTS
} ParentStatus;

typedef enum {
  TREE_STATUS_UNDETERMINED,
  IS_TREE,
  IS_NOT_A_TREE
} TreeStatus;

typedef struct {
  GtGenomeNodeTraverseFunc func;
  void *data;
} Traverse_children_info;

typedef struct {
  GtArray *exon_features,
          *cds_features;
} SaveExonAndCDSInfo;

static void gt_feature_node_free(GtGenomeNode *gn)
{
  GtFeatureNode *fn = gt_feature_node_cast(gn);
  gt_str_delete(fn->seqid);
  gt_str_delete(fn->source);
  gt_tag_value_map_delete(fn->attributes);
  gt_dlist_delete(fn->children);
}

const char* gt_feature_node_get_attribute(GtFeatureNode *fn,
                                          const char *attr_name)
{
  if (!fn->attributes)
    return NULL;
  return gt_tag_value_map_get(fn->attributes, attr_name);
}

static void store_attribute(const char *attr_name,
                            GT_UNUSED const char *attr_value, void *data)
{
  GtStrArray *list = data;
  gt_assert(attr_name && attr_value && data);
  gt_str_array_add_cstr(list, attr_name);
}

GtStrArray* gt_feature_node_get_attribute_list(GtFeatureNode *fn)
{
  GtStrArray *list = gt_str_array_new();
  if (fn->attributes)
    gt_tag_value_map_foreach(fn->attributes, store_attribute, list);
  return list;
}

static GtStr* gt_feature_node_get_seqid(GtGenomeNode *gn)
{
  GtFeatureNode *fn = gt_feature_node_cast(gn);
  return fn->seqid;
}

static GtRange gt_feature_node_get_range(GtGenomeNode *gn)
{
  GtFeatureNode *fn = gt_feature_node_cast(gn);
  return fn->range;
}

static void gt_feature_node_set_range(GtGenomeNode *gn, const GtRange *range)
{
  GtFeatureNode *fn = gt_feature_node_cast(gn);
  fn->range = *range;
}

static void gt_feature_node_change_seqid(GtGenomeNode *gn, GtStr *seqid)
{
  GtFeatureNode *fn = gt_feature_node_cast(gn);
  gt_assert(fn && seqid);
  gt_str_delete(fn->seqid);
  fn->seqid = gt_str_ref(seqid);
}

void gt_feature_node_set_source(GtFeatureNode *fn, GtStr *source)
{
  gt_assert(fn && source && !fn->source);
  fn->source = gt_str_ref(source);
}

void gt_feature_node_set_phase(GtFeatureNode *fn, GtPhase phase)
{
  gt_assert(fn);
  fn->bit_field &= ~(PHASE_MASK << PHASE_OFFSET);
  fn->bit_field |= phase << PHASE_OFFSET;
}

static int gt_feature_node_accept(GtGenomeNode *gn, GtNodeVisitor *gv,
                                  GtError *err)
{
  GtFeatureNode *fn;
  gt_error_check(err);
  fn = gt_feature_node_cast(gn);
  return gt_node_visitor_visit_feature_node(gv, fn, err);
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

static void set_transcriptfeaturetype(GtFeatureNode *fn,
                                      GtTranscriptFeatureType tft)
{
  gt_assert(fn);
  fn->bit_field &= ~(TRANSCRIPT_FEATURE_TYPE_MASK <<
                     TRANSCRIPT_FEATURE_TYPE_OFFSET);
  fn->bit_field |= tft << TRANSCRIPT_FEATURE_TYPE_OFFSET;
}

static void set_tree_status(unsigned int *bit_field, TreeStatus tree_status)
{
  *bit_field &= ~(TREE_STATUS_MASK << TREE_STATUS_OFFSET);
  *bit_field |= tree_status << TREE_STATUS_OFFSET;
}

GtGenomeNode* gt_feature_node_new(GtStr *seqid, const char *type,
                                  unsigned long start, unsigned long end,
                                  GtStrand strand)
{
  GtGenomeNode *gn;
  GtFeatureNode *fn;
  gt_assert(seqid && type);
  gt_assert(start <= end);
  gn = gt_genome_node_create(gt_feature_node_class());
  fn = gt_feature_node_cast(gn);
  fn->seqid       = gt_str_ref(seqid);
  fn->source      = NULL;
  fn->type        = gt_symbol(type);
  fn->score       = UNDEF_FLOAT;
  fn->range.start = start;
  fn->range.end   = end;
  fn->attributes  = NULL;
  fn->bit_field   = 0;
  fn->bit_field |= strand << STRAND_OFFSET;
  fn->children    = NULL; /* the children list is create on demand */
  gt_feature_node_set_phase(fn, GT_PHASE_UNDEFINED);
  set_transcriptfeaturetype(fn, TRANSCRIPT_FEATURE_TYPE_UNDETERMINED);
  set_tree_status(&fn->bit_field, IS_TREE);
  fn->representative = NULL;
  return gn;
}

GtGenomeNode* gt_feature_node_new_pseudo(GtFeatureNode *fn)
{
  GtFeatureNode *pf;
  GtGenomeNode *pn;
  GtRange range;
  gt_assert(fn);
  range = gt_feature_node_get_range((GtGenomeNode*) fn),
  pn = gt_feature_node_new(gt_feature_node_get_seqid((GtGenomeNode*) fn),
                           gt_feature_node_get_type(fn), range.start,
                           range.end, gt_feature_node_get_strand(fn));
  pf = gt_feature_node_cast(pn);
  pf->type = NULL; /* pseudo features do not have a type */
  gt_feature_node_set_source(pf, fn->source);
  pf->bit_field |= 1 << PSEUDO_FEATURE_OFFSET;
  return pn;
}

GtGenomeNode* gt_feature_node_new_standard_gene(void)
{
  GtGenomeNode *fn, *child, *grandch;
  GtStr *seqid;
  seqid = gt_str_new_cstr("ctg123");

  /* gene */
  fn = gt_feature_node_new(seqid, gft_gene, 1000, 9000, GT_STRAND_FORWARD);

  /* TF binding site */
  child = gt_feature_node_new(seqid, gft_TF_binding_site, 1000, 1012,
                                GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) fn, (GtFeatureNode*) child);

  /* first mRNA */
  child = gt_feature_node_new(seqid, gft_mRNA, 1050, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) fn, (GtFeatureNode*) child);
  grandch = gt_feature_node_new(seqid, gft_exon, 1050, 1500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grandch);
  grandch = gt_feature_node_new(seqid, gft_exon, 3000, 3902, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grandch);
  grandch = gt_feature_node_new(seqid, gft_exon, 5000, 5500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grandch);
  grandch = gt_feature_node_new(seqid, gft_exon, 7000, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grandch);

  /* second mRNA */
  child = gt_feature_node_new(seqid, gft_mRNA, 1050, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) fn, (GtFeatureNode*) child);
  grandch = gt_feature_node_new(seqid, gft_exon, 1050, 1500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grandch);
  grandch = gt_feature_node_new(seqid, gft_exon, 5000, 5500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grandch);
  grandch = gt_feature_node_new(seqid, gft_exon, 7000, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grandch);

  /* third mRNA */
  child = gt_feature_node_new(seqid, gft_mRNA, 1300, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) fn, (GtFeatureNode*) child);
  grandch = gt_feature_node_new(seqid, gft_exon, 1300, 1500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grandch);
  grandch = gt_feature_node_new(seqid, gft_exon, 3000, 3902, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grandch);
  grandch = gt_feature_node_new(seqid, gft_exon, 5000, 5500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grandch);
  grandch = gt_feature_node_new(seqid, gft_exon, 7000, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grandch);

  gt_str_delete(seqid);
  return fn;
}

const char* gt_feature_node_get_source(GtFeatureNode *fn)
{
  gt_assert(fn);
  return fn->source ? gt_str_get(fn->source) : ".";
}

const char* gt_feature_node_get_type(GtFeatureNode *fn)
{
  gt_assert(fn);
  return fn->type;
}

bool gt_feature_node_has_type(GtFeatureNode *fn, const char *type)
{
  gt_assert(fn && type);
  return gt_strcmp(fn->type, type) ? false : true;
}

bool gt_feature_node_score_is_defined(const GtFeatureNode *fn)
{
  gt_assert(fn);
  if ((fn->bit_field >> SCORE_IS_DEFINED_OFFSET) & SCORE_IS_DEFINED_MASK)
    return true;
  return false;
}

bool gt_feature_node_is_multi(const GtFeatureNode *fn)
{
  gt_assert(fn);
  if ((fn->bit_field >> MULTI_FEATURE_OFFSET) & MULTI_FEATURE_MASK) {
    gt_assert(!gt_feature_node_is_pseudo(fn));
    return true;
  }
  return false;
}

bool gt_feature_node_is_pseudo(const GtFeatureNode *fn)
{
  gt_assert(fn);
  if ((fn->bit_field >> PSEUDO_FEATURE_OFFSET) & PSEUDO_FEATURE_MASK) {
    gt_assert(!gt_feature_node_is_multi(fn));
    return true;
  }
  return false;
}

static void gt_feature_node_set_multi(GtFeatureNode *fn)
{
  gt_assert(fn && !gt_feature_node_is_multi(fn));
  fn->bit_field |= 1 << MULTI_FEATURE_OFFSET;
}

void gt_feature_node_make_multi_representative(GtFeatureNode *fn)
{
  gt_assert(fn && !gt_feature_node_is_multi(fn));
  gt_feature_node_set_multi(fn);
}

void gt_feature_node_set_multi_representative(GtFeatureNode *fn,
                                             GtFeatureNode *rep)
{
  gt_assert(fn && !gt_feature_node_is_multi(fn));
  gt_assert(rep && gt_feature_node_is_multi(rep));
  gt_feature_node_set_multi(fn);
  fn->representative = rep;
}

GtFeatureNode* gt_feature_node_get_multi_representative(GtFeatureNode
                                                             *fn)
{
  gt_assert(fn && gt_feature_node_is_multi(fn) &&
         !gt_feature_node_is_pseudo(fn));
  if (fn->representative) {
    gt_assert(gt_feature_node_is_multi(fn->representative));
    gt_assert(gt_feature_node_get_multi_representative(fn->representative) ==
           fn->representative);
    return fn->representative;
  }
  return fn; /* is itself the representative */
}

float gt_feature_node_get_score(const GtFeatureNode *fn)
{
  gt_assert(fn);
  gt_assert(gt_feature_node_score_is_defined(fn));
  return fn->score;
}

GtStrand gt_feature_node_get_strand(GtFeatureNode *fn)
{
  gt_assert(fn);
  return (fn->bit_field >> STRAND_OFFSET) & STRAND_MASK;
}

void gt_feature_node_set_strand(GtFeatureNode *fn, GtStrand strand)
{
  gt_assert(fn);
  fn->bit_field &= ~(STRAND_MASK << STRAND_OFFSET);
  fn->bit_field |= strand << STRAND_OFFSET;
}

GtPhase gt_feature_node_get_phase(GtFeatureNode *fn)
{
  gt_assert(fn);
  return (fn->bit_field >> PHASE_OFFSET) & PHASE_MASK;
}

static int feature_node_save_exon(GtGenomeNode *gn, void *data,
                                  GT_UNUSED GtError *err)
{
  GtFeatureNode *fn;
  GtArray *exon_features = (GtArray*) data;
  gt_error_check(err);
  fn = (GtFeatureNode*) gn;
  gt_assert(fn && exon_features);
  if (gt_feature_node_has_type(fn, gft_exon)) {
    gt_array_add(exon_features, fn);
  }
  return 0;
}

void gt_feature_node_get_exons(GtFeatureNode *fn, GtArray *exon_features)
{
  int had_err;
  gt_assert(fn && exon_features && !gt_array_size(exon_features));
  had_err = gt_genome_node_traverse_children((GtGenomeNode*) fn, exon_features,
                                             feature_node_save_exon, false,
                                             NULL);
  gt_assert(!had_err); /* feature_node_save_exon() is sane */
}

static int save_exons_and_cds(GtGenomeNode *gn, void *data,
                              GT_UNUSED GtError *err)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  GtFeatureNode *fn;
  gt_error_check(err);
  fn = (GtFeatureNode*) gn;
  gt_assert(fn && info);
  if (gt_feature_node_has_type(fn, gft_exon))
    gt_array_add(info->exon_features, fn);
  else if (gt_feature_node_has_type(fn, gft_CDS))
    gt_array_add(info->cds_features, fn);
  return 0;
}

static void set_transcript_types(GtArray *features)
{
  GtFeatureNode *fn;
  unsigned long i;
  gt_assert(features);
  if (gt_array_size(features)) {
    if (gt_array_size(features) == 1) {
      fn = *(GtFeatureNode**) gt_array_get(features, 0);
      set_transcriptfeaturetype(fn, TRANSCRIPT_FEATURE_TYPE_SINGLE);
    }
    else {
      fn = *(GtFeatureNode**) gt_array_get(features, 0);
      set_transcriptfeaturetype(fn, TRANSCRIPT_FEATURE_TYPE_INITIAL);
      for (i = 1; i < gt_array_size(features) - 1; i++) {
        fn = *(GtFeatureNode**) gt_array_get(features, i);
        set_transcriptfeaturetype(fn, TRANSCRIPT_FEATURE_TYPE_INTERNAL);
      }
      fn = *(GtFeatureNode**)
           gt_array_get(features, gt_array_size(features) - 1);
      set_transcriptfeaturetype(fn, TRANSCRIPT_FEATURE_TYPE_TERMINAL);
    }
  }
}

static int determine_transcripttypes(GtGenomeNode *gn, void *data,
                                     GT_UNUSED GtError *err)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  int had_err;
  gt_error_check(err);
  gt_assert(gn && info);
  /* reset exon_features and cds_features */
  gt_array_reset(info->exon_features);
  gt_array_reset(info->cds_features);
  /* collect all direct children exons */
  had_err = gt_genome_node_traverse_direct_children(gn, info,
                                                    save_exons_and_cds, NULL);
  gt_assert(!had_err); /* cannot happen, because save_exons_and_cds() is sane */
  /* set transcript feature type, if necessary */
  set_transcript_types(info->exon_features);
  set_transcript_types(info->cds_features);
  return 0;
}

void gt_feature_node_determine_transcripttypes(GtFeatureNode *fn)
{
  SaveExonAndCDSInfo info;
  int had_err;
  gt_assert(fn);
  info.exon_features = gt_array_new(sizeof (GtFeatureNode*));
  info.cds_features = gt_array_new(sizeof (GtFeatureNode*));
  had_err = gt_genome_node_traverse_children((GtGenomeNode*) fn, &info,
                                          determine_transcripttypes, false,
                                          NULL);
  gt_assert(!had_err); /* cannot happen, because determine_transcripttypes() is
                       sane */
  gt_array_delete(info.exon_features);
  gt_array_delete(info.cds_features);
}

GtTranscriptFeatureType gt_feature_node_get_transcriptfeaturetype(GtFeatureNode
                                                                *fn)
{
  gt_assert(fn);
  return (fn->bit_field >> TRANSCRIPT_FEATURE_TYPE_OFFSET) &
         TRANSCRIPT_FEATURE_TYPE_MASK;
}

void gt_feature_node_set_end(GtFeatureNode *fn, unsigned long end)
{
  gt_assert(fn && fn->range.start <= end);
  fn->range.end = end;
}

void gt_feature_node_set_score(GtFeatureNode *fn, float score)
{
  gt_assert(fn);
  fn->bit_field |= 1 << SCORE_IS_DEFINED_OFFSET;
  fn->score = score;
}

void gt_feature_node_unset_score(GtFeatureNode *fn)
{
  gt_assert(fn);
  fn->bit_field &= ~(1 << SCORE_IS_DEFINED_OFFSET);
  fn->score = UNDEF_FLOAT;
}

void gt_feature_node_add_attribute(GtFeatureNode *fn,
                                   const char *attr_name,
                                   const char *attr_value)
{
  gt_assert(fn && attr_name && attr_value);
  gt_assert(strlen(attr_name)); /* attribute name cannot be empty */
  gt_assert(strlen(attr_value)); /* attribute value cannot be empty */
  if (!fn->attributes)
    fn->attributes = gt_tag_value_map_new(attr_name, attr_value);
  else
    gt_tag_value_map_add(&fn->attributes, attr_name, attr_value);
}

void gt_feature_node_foreach_attribute(GtFeatureNode *fn,
                                      AttributeIterFunc iterfunc, void *data)
{
  gt_assert(fn && iterfunc);
  if (fn->attributes) {
    gt_tag_value_map_foreach(fn->attributes,
                             (GtTagValueMapIteratorFunc) iterfunc,
                             data);
  }
}

static bool feature_node_has_gft(const GtFeatureNode *fn, const char **fnts)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *child;
  bool has_gft = false;
  gt_assert(fn && fnts && fnts[0] != NULL);
  fni = gt_feature_node_iterator_new(fn);
  while ((child = gt_feature_node_iterator_next(fni))) {
    unsigned long i = 0;
    while (fnts[i] != NULL) {
      if (gt_feature_node_has_type((GtFeatureNode*) child, fnts[i])) {
        has_gft = true;
        break;
      }
      i++;
    }
    if (has_gft)
      break;
  }
  gt_feature_node_iterator_delete(fni);
  return has_gft;
}

bool gt_feature_node_has_CDS(const GtFeatureNode *fn)
{
  static const char *gfts[] = { gft_CDS, NULL };
  return feature_node_has_gft(fn, gfts);
}

bool gt_feature_node_has_splice_site(const GtFeatureNode *fn)
{
  static const char *gfts[] = { gft_five_prime_splice_site,
                                gft_three_prime_splice_site, NULL };
  return feature_node_has_gft(fn, gfts);
}

double gt_feature_node_average_splice_site_prob(const GtFeatureNode *fn)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *child;
  unsigned long num_of_splice_sites = 0;
  double averagessp = 0.0;
  gt_assert(fn);
  fni = gt_feature_node_iterator_new(fn);
  while ((child = gt_feature_node_iterator_next(fni))) {
    if (gt_feature_node_has_type(child, gft_five_prime_splice_site) ||
        gt_feature_node_has_type(child, gft_three_prime_splice_site)) {
      averagessp += gt_feature_node_get_score((GtFeatureNode*) child);
      num_of_splice_sites++;
    }
  }
  gt_feature_node_iterator_delete(fni);
  if (num_of_splice_sites)
    averagessp /= num_of_splice_sites;
  return averagessp;
}

bool gt_genome_features_are_similar(GtFeatureNode *fn_a,
                                    GtFeatureNode *fn_b)
{
  GtRange range_a, range_b;
  gt_assert(fn_a && fn_b);
  range_a = gt_feature_node_get_range((GtGenomeNode*) fn_a);
  range_b = gt_feature_node_get_range((GtGenomeNode*) fn_b);
  if (!gt_str_cmp(gt_genome_node_get_seqid((GtGenomeNode*) fn_a),
                  gt_genome_node_get_seqid((GtGenomeNode*) fn_b)) &&
      (gt_feature_node_get_type(fn_a) == gt_feature_node_get_type(fn_b)) &&
      (!gt_range_compare(&range_a, &range_b)) &&
      (gt_feature_node_get_strand(fn_a) ==
       gt_feature_node_get_strand(fn_b)) &&
      (gt_feature_node_get_phase(fn_a) ==
       gt_feature_node_get_phase(fn_b))) {
    return true;
  }
  return false;
}

int gt_feature_node_unit_test(GtError *err)
{
  GtGenomeNode *fn;
  GtStr *seqid;
  int had_err = 0;

  gt_error_check(err);

  seqid = gt_str_new_cstr("seqid");
  fn = gt_feature_node_new(seqid, gft_gene, 1, 1000, GT_STRAND_FORWARD);

  ensure(had_err, !gt_feature_node_score_is_defined((GtFeatureNode*) fn));

  gt_genome_node_delete(fn);
  gt_str_delete(seqid);

  return had_err;
}

static TreeStatus get_tree_status(unsigned int bit_field)
{
  return (bit_field >> TREE_STATUS_OFFSET) & TREE_STATUS_MASK;
}

static void set_parent_status(unsigned int *bit_field,
                              ParentStatus parent_status)
{
  *bit_field &= ~(PARENT_STATUS_MASK << PARENT_STATUS_OFFSET);
  *bit_field |= parent_status << PARENT_STATUS_OFFSET;
}

static ParentStatus get_parent_status(unsigned int bit_field)
{
  return (bit_field >> PARENT_STATUS_OFFSET) & PARENT_STATUS_MASK;
}

static bool multiple_parents(unsigned int bit_field)
{
  if (get_parent_status(bit_field) == MULTIPLE_PARENTS)
    return true;
  return false;
}

static void add_parent(unsigned int *bit_field)
{
  switch (get_parent_status(*bit_field)) {
    case NO_PARENT:
      set_parent_status(bit_field, ONE_PARENT);
      break;
    case ONE_PARENT:
      set_parent_status(bit_field, MULTIPLE_PARENTS);
      break;
    case MULTIPLE_PARENTS:
      break;
  }
}

/* XXX: static? */
int gt_genome_node_traverse_children_generic(GtGenomeNode *genome_node,
                                          void *data,
                                          GtGenomeNodeTraverseFunc traverse,
                                          bool traverse_only_once,
                                          bool depth_first, bool with_pseudo,
                                          GtError *err)
{
  GtArray *node_stack = NULL, *list_of_children;
  GtQueue *node_queue = NULL;
  GtGenomeNode *gn, *child_feature;
  GtFeatureNode *feature_node, *fn, *fn_ref;
  GtDlistelem *dlistelem;
  unsigned long i;
  GtHashtable *traversed_nodes = NULL;
  bool has_node_with_multiple_parents = false;
  int had_err = 0;

  if (!genome_node)
    return 0;

  /* XXX */
  feature_node = gt_feature_node_cast(genome_node);

  /* create additional reference to <genome_node> (necessary if genome_node is
     freed by <traverse>) */
  fn_ref = gt_feature_node_nonrec_ref(feature_node);

  if (depth_first) {
    node_stack = gt_array_new(sizeof (GtGenomeNode*));
    if (!with_pseudo && gt_feature_node_try_cast(genome_node) &&
        gt_feature_node_is_pseudo((GtFeatureNode*) genome_node)) {
      /* add the children backwards to traverse in order */
      for (dlistelem = gt_dlist_last(feature_node->children);
           dlistelem != NULL;
           dlistelem = gt_dlistelem_previous(dlistelem)) {
        child_feature = (GtGenomeNode*) gt_dlistelem_get_data(dlistelem);
        gt_array_add(node_stack, child_feature);
      }
    }
    else
      gt_array_add(node_stack, feature_node);
    gt_assert(gt_array_size(node_stack));
  }
  else {
    node_queue = gt_queue_new();
    if (!with_pseudo && gt_feature_node_is_pseudo(feature_node)) {
      for (dlistelem = gt_dlist_first(feature_node->children);
           dlistelem != NULL;
           dlistelem = gt_dlistelem_next(dlistelem)) {
        child_feature = (GtGenomeNode*) gt_dlistelem_get_data(dlistelem);
        gt_queue_add(node_queue, child_feature);
      }
    }
    else
      gt_queue_add(node_queue, feature_node);
    gt_assert(gt_queue_size(node_queue));
  }
  list_of_children = gt_array_new(sizeof (GtGenomeNode*));

  if (traverse_only_once)
  {
    static const HashElemInfo node_hashtype
      = { gt_ht_ptr_elem_hash, { NULL }, sizeof (GtGenomeNode *),
          gt_ht_ptr_elem_cmp, NULL, NULL };
    traversed_nodes = gt_hashtable_new(node_hashtype);
  }

  while ((depth_first ? gt_array_size(node_stack)
                      : gt_queue_size(node_queue))) {
    if (depth_first)
      gn = *(GtGenomeNode**) gt_array_pop(node_stack);
    else
      gn = gt_queue_get(node_queue);
    gt_array_reset(list_of_children);
    /* XXX */
    fn = gt_feature_node_cast(gn);
    if (fn->children) {
      /* a backup of the children array is necessary if traverse() frees the
         node */
      for (dlistelem = gt_dlist_first(fn->children); dlistelem != NULL;
           dlistelem = gt_dlistelem_next(dlistelem)) {
        child_feature = (GtGenomeNode*) gt_dlistelem_get_data(dlistelem);
        gt_array_add(list_of_children, child_feature);
      }
    }
    /* store the implications of <gn> to the tree status of <feature_node> */
    if (multiple_parents(fn->bit_field))
      has_node_with_multiple_parents = true;
    /* call traverse function */
    if (traverse) {
      had_err = traverse(gn, data, err);
      if (had_err)
        break;
    }
    for (i = 0; i < gt_array_size(list_of_children); i++) {
      if (depth_first) {
        /* we go backwards to traverse in order */
        child_feature = *(GtGenomeNode**) gt_array_get(list_of_children,
                                       gt_array_size(list_of_children) - i - 1);
      }
      else {
        child_feature = *(GtGenomeNode**) gt_array_get(list_of_children, i);
      }
      if (!traverse_only_once ||
          !gt_hashtable_get(traversed_nodes, &child_feature)) {
        /* feature has not been traversed or has to be traversed multiple
           times */
        if (depth_first)
          gt_array_add(node_stack, child_feature);
        else
          gt_queue_add(node_queue, child_feature);
        if (traverse_only_once)
          gt_hashtable_add(traversed_nodes, &child_feature);
      }
    }
  }

  /* save the tree status of the genome node */
  if (!had_err) {
    if (has_node_with_multiple_parents) {
      set_tree_status(&fn_ref->bit_field, IS_NOT_A_TREE);
      gt_assert(get_tree_status(fn_ref->bit_field) == IS_NOT_A_TREE);
    }
    else {
      set_tree_status(&fn_ref->bit_field, IS_TREE);
      gt_assert(get_tree_status(fn_ref->bit_field) ==IS_TREE);
    }
  }

  /* free */
  gt_feature_node_nonrec_delete(fn_ref);
  if (traverse_only_once)
    gt_hashtable_delete(traversed_nodes);
  gt_array_delete(list_of_children);
  gt_array_delete(node_stack);
  gt_queue_delete(node_queue);

  return had_err;
}

static int gt_genome_node_traverse_children_with_pseudo(GtGenomeNode
                                                        *genome_node,
                                                        void *data,
                                                       GtGenomeNodeTraverseFunc
                                                        traverse,
                                                        bool traverse_only_once,
                                                        GtError *err)
{
  return gt_genome_node_traverse_children_generic(genome_node, data, traverse,
                                               traverse_only_once, true, true,
                                               err);
}

static int increase_reference_count(GtGenomeNode *gn, GT_UNUSED void *data,
                                    GT_UNUSED GtError *err)
{
  gt_error_check(err);
  gt_assert(gn);
  gn->reference_count++;
  return 0;
}

GtGenomeNode* gt_genome_node_ref(GtGenomeNode *gn)
{
  GtFeatureNode *fn;
  int had_err;
  gt_assert(gn);
  if ((fn = gt_feature_node_try_cast(gn))) {
    had_err = gt_genome_node_traverse_children_with_pseudo(gn, NULL,
                                                      increase_reference_count,
                                                           true, NULL);
    gt_assert(!had_err); /* cannot happen, increase_reference_count() is sane */
  }
  else
    gn->reference_count++;
  return gn;
}

int gt_genome_node_traverse_children(GtGenomeNode *genome_node, void *data,
                                  GtGenomeNodeTraverseFunc traverse,
                                  bool traverse_only_once, GtError *err)
{
  return gt_genome_node_traverse_children_generic(genome_node, data, traverse,
                                               traverse_only_once, true, false,
                                               err);
}

int gt_genome_node_traverse_children_breadth(GtGenomeNode *genome_node,
                                             void *data,
                                             GtGenomeNodeTraverseFunc traverse,
                                             bool traverse_only_once,
                                             GtError *err)
{
  return gt_genome_node_traverse_children_generic(genome_node, data, traverse,
                                               traverse_only_once, false, false,
                                               err);
}

int gt_genome_node_traverse_direct_children(GtGenomeNode *gn,
                                         void *traverse_func_data,
                                         GtGenomeNodeTraverseFunc traverse,
                                         GtError *err)
{
  GtFeatureNode *fn;
  GtDlistelem *dlistelem;
  int had_err = 0;
  gt_error_check(err);
  if (!gn || !traverse)
    return 0;
  /* XXX */
  fn = gt_feature_node_cast(gn);
  if (fn->children) {
    for (dlistelem = gt_dlist_first(fn->children); dlistelem != NULL;
         dlistelem = gt_dlistelem_next(dlistelem)) {
      had_err = traverse((GtGenomeNode*) gt_dlistelem_get_data(dlistelem),
                          traverse_func_data, err);
      if (had_err)
        break;
    }
  }
  return had_err;
}

unsigned long gt_genome_node_number_of_children(const GtGenomeNode *gn)
{
  GtFeatureNode *fn;
  gt_assert(gn);
  fn = gt_feature_node_cast((GtGenomeNode*) gn); /* XXX */
  return gt_dlist_size(fn->children);
}

void gt_feature_node_add_child(GtFeatureNode *parent, GtFeatureNode *child)
{
  gt_assert(parent && child);
  /* <parent> and <child> have the same seqid */
  gt_assert(!gt_str_cmp(gt_genome_node_get_seqid((GtGenomeNode*) parent),
                        gt_genome_node_get_seqid((GtGenomeNode*) child)));
  /* pseudo-features have to be top-level */
  gt_assert(!gt_feature_node_is_pseudo((GtFeatureNode*) child));
  /* create children list on demand */
  if (!parent->children)
    parent->children = gt_dlist_new((GtCompare) gt_genome_node_cmp);
  gt_dlist_add(parent->children, child); /* XXX: check for circles */
  /* update tree status of <parent> */
  set_tree_status(&parent->bit_field, TREE_STATUS_UNDETERMINED);
  /* update parent info of <child> */
  add_parent(&child->bit_field);
}

static int remove_leaf(GtGenomeNode *node, void *data, GT_UNUSED GtError *err)
{
  GtFeatureNode *node_feature;
  GtDlistelem *dlistelem;
  GtGenomeNode *child, *leaf = (GtGenomeNode*) data;
  gt_error_check(err);
  node_feature = gt_feature_node_cast(node); /* XXX */
  if (node != leaf && node_feature->children) {
    for (dlistelem = gt_dlist_first(node_feature->children); dlistelem != NULL;
         dlistelem = gt_dlistelem_next(dlistelem)) {
      child = (GtGenomeNode*) gt_dlistelem_get_data(dlistelem);
      if (child == leaf) {
        gt_dlist_remove(node_feature->children, dlistelem);
        break;
      }
    }
  }
  return 0;
}

void gt_genome_node_remove_leaf(GtGenomeNode *tree, GtGenomeNode *leafn)
{
  int had_err;
  gt_assert(tree && leafn);
  gt_assert(!gt_genome_node_number_of_children(leafn));
  had_err = gt_genome_node_traverse_children(tree, leafn, remove_leaf, true,
                                             NULL);
  gt_assert(!had_err); /* cannot happen, remove_leaf() is sane */
}

void gt_genome_node_mark(GtGenomeNode *gn)
{
  GtFeatureNode *fn;
  gt_assert(gn);
  fn = gt_feature_node_cast(gn); /* XXX */
  fn->bit_field |= 1;
}

bool gt_genome_node_is_marked(const GtGenomeNode *gn)
{
  GtFeatureNode *fn;
  gt_assert(gn);
  fn = gt_feature_node_cast((GtGenomeNode*) gn); /* XXX */
  return fn->bit_field & 1 ? true : false;
}

static int check_marked_status(GtGenomeNode *gn, void *data,
                               GT_UNUSED GtError *err)
{
  bool *marked = data;
  if (gt_genome_node_is_marked(gn))
    *marked = true;
  return 0;
}

bool gt_genome_node_contains_marked(GtGenomeNode *gn)
{
  bool contains_marked = false;
  int rval;
  gt_assert(gn);
  rval = gt_genome_node_traverse_children(gn, &contains_marked,
                                       check_marked_status, true, NULL);
  gt_assert(!rval); /* check_marked_status() is sane */
  return contains_marked;
}

bool gt_genome_node_has_children(GtGenomeNode *gn)
{
  GtFeatureNode *fn;
  gt_assert(gn);
  fn = gt_feature_node_cast((GtGenomeNode*) gn); /* XXX */
  if (!fn->children || gt_dlist_size(fn->children) == 0)
    return false;
  return true;
}

bool gt_genome_node_direct_children_do_not_overlap_generic(GtGenomeNode
                                                           *parent,
                                                           GtGenomeNode
                                                           *child)
{
  GtFeatureNode *parent_node;
  GtArray *children_ranges;
  GtDlistelem *dlistelem;
  GtFeatureNode *fn = NULL, *child_fn;
  GtRange range;
  bool rval;

  gt_assert(parent);

  if (child)
    fn = gt_feature_node_try_cast(child);

  parent_node = gt_feature_node_cast(parent); /* XXX */

  if (!parent_node->children)
    return true;

  /* get children ranges */
  children_ranges = gt_array_new(sizeof (GtRange));
  gt_assert(parent_node->children);
  for (dlistelem = gt_dlist_first(parent_node->children); dlistelem != NULL;
       dlistelem = gt_dlistelem_next(dlistelem)) {
    if (!fn ||
        ((child_fn =
            gt_feature_node_try_cast(gt_dlistelem_get_data(dlistelem))) &&
         gt_feature_node_get_type(fn) ==
         gt_feature_node_get_type(child_fn))) {
      range = gt_genome_node_get_range((GtGenomeNode*)
                                    gt_dlistelem_get_data(dlistelem));
      gt_array_add(children_ranges, range);
    }
  }

  gt_ranges_sort(children_ranges);
  gt_assert(gt_ranges_are_sorted(children_ranges));
  rval = gt_ranges_do_not_overlap(children_ranges);

  gt_array_delete(children_ranges);

  return rval;
}

bool gt_genome_node_direct_children_do_not_overlap(GtGenomeNode *gn)
{
  return gt_genome_node_direct_children_do_not_overlap_generic(gn, NULL);
}

bool gt_genome_node_direct_children_do_not_overlap_st(GtGenomeNode *parent,
                                                   GtGenomeNode *child)
{
  return gt_genome_node_direct_children_do_not_overlap_generic(parent, child);
}

bool gt_genome_node_is_tree(GtGenomeNode *gn)
{
  GtFeatureNode *fn;
  bool status = false;
  fn = gt_feature_node_cast(gn); /* XXX */
  gt_assert(gn);
  switch (get_tree_status(fn->bit_field)) {
    case IS_TREE:
      status = true;
      break;
    case IS_NOT_A_TREE:
      status = false;
      break;
    case TREE_STATUS_UNDETERMINED:
      /* not implemented, the tree status must have been determined by a
         previous gt_genome_node_traverse_children() invocation */
    default: gt_assert(0);
  }
  return status;
}

bool gt_genome_node_overlaps_nodes(GtGenomeNode *gn, GtArray *nodes)
{
  return gt_genome_node_overlaps_nodes_mark(gn, nodes, NULL);
}

bool gt_genome_node_overlaps_nodes_mark(GtGenomeNode *gn, GtArray *nodes,
                                             GtBittab *b)
{
  unsigned long i;
  GtGenomeNode *node;
  GtRange gn_range, node_range;
  bool rval = false;
#ifndef NDEBUG
  GtStr *gn_id;
  gt_assert(gn && nodes);
  gt_assert(!b || gt_bittab_size(b) == gt_array_size(nodes));
  gn_id = gt_genome_node_get_idstr(gn);
#endif
  gn_range = gt_genome_node_get_range(gn);

  for (i = 0; i < gt_array_size(nodes); i++) {
    node = *(GtGenomeNode**) gt_array_get(nodes, i);
    gt_assert(!gt_str_cmp(gn_id, gt_genome_node_get_idstr(node)));
    node_range = gt_genome_node_get_range(node);
    if (gt_range_overlap(&gn_range, &node_range)) {
      rval = true;
      if (b)
        gt_bittab_set_bit(b, i);
      else
        break;
    }
  }
  return rval;
}

static int free_genome_node(GtGenomeNode *gn, GT_UNUSED void *data,
                            GT_UNUSED GtError *err)
{
  gt_feature_node_nonrec_delete((GtFeatureNode*) gn);
  return 0;
}

void gt_genome_node_delete(GtGenomeNode *gn)
{
  GtFeatureNode *fn;
  int had_err;
  if (!gn) return;
  /* XXX */
  if ((fn = gt_feature_node_try_cast(gn))) {
    had_err = gt_genome_node_traverse_children_with_pseudo(gn, NULL,
                                                         free_genome_node, true,
                                                           NULL);
    gt_assert(!had_err); /* cannot happen, free_genome_node() is sane */
  }
  else {
    if (gn->reference_count) {
      gn->reference_count--;
      return;
    }
    gt_assert(gn->c_class);
    if (gn->c_class->free)
      gn->c_class->free(gn);
    gt_str_delete(gn->filename);
    if (gn->userdata)
      gt_hashmap_delete(gn->userdata);
    gt_free(gn);
  }
}

GtFeatureNode* gt_feature_node_nonrec_ref(GtFeatureNode *fn)
{
  gt_assert(fn);
  fn->parent_instance.reference_count++;
  return fn;
}

void gt_feature_node_nonrec_delete(GtFeatureNode *fn)
{
  if (!fn) return;
  if (fn->parent_instance.reference_count) {
    fn->parent_instance.reference_count--;
    return;
  }
  gt_assert(fn->parent_instance.c_class);
  if (fn->parent_instance.c_class->free)
    fn->parent_instance.c_class->free((GtGenomeNode*) fn);
  gt_str_delete(fn->parent_instance.filename);
  if (fn->parent_instance.userdata)
      gt_hashmap_delete(fn->parent_instance.userdata);
  gt_free(fn);
}
