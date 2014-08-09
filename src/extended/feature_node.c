/*
  Copyright (c) 2006-2013 Gordon Gremme <gordon@gremme.org>
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
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/ensure.h"
#include "core/hashtable.h"
#include "core/ma.h"
#include "core/queue_api.h"
#include "core/strcmp_api.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_rep.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/genome_node_rep.h"

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
#define DFS_STATUS_OFFSET               16
#define DFS_STATUS_MASK                 0x3

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

typedef enum {
  DFS_WHITE,
  DFS_GRAY,
  DFS_BLACK
} DFSStatus;

typedef struct {
  GtArray *exon_features,
          *cds_features;
} SaveExonAndCDSInfo;

typedef struct {
  const char *type;
  GtUword number;
} GtTypeTraverseInfo;

static void feature_node_free(GtGenomeNode *gn)
{
  GtFeatureNode *fn = gt_feature_node_cast(gn);
  gt_str_delete(fn->seqid);
  gt_str_delete(fn->source);
  gt_tag_value_map_delete(fn->attributes);
  if (fn->children) {
    GtDlistelem *dlistelem;
    for (dlistelem = gt_dlist_first(fn->children);
         dlistelem != NULL;
         dlistelem = gt_dlistelem_next(dlistelem)) {
      gt_genome_node_delete(gt_dlistelem_get_data(dlistelem));
    }
  }
  gt_dlist_delete(fn->children);
  if (fn->observer && fn->observer->deleted)
    fn->observer->deleted(fn, fn->observer->data);
  if (fn->observer)
    gt_feature_node_observer_delete(fn->observer);
}

const char* gt_feature_node_get_attribute(const GtFeatureNode *fn,
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

GtStrArray* gt_feature_node_get_attribute_list(const GtFeatureNode *fn)
{
  GtStrArray *list = gt_str_array_new();
  if (fn->attributes)
    gt_tag_value_map_foreach(fn->attributes, store_attribute, list);
  return list;
}

static GtStr* feature_node_get_seqid(GtGenomeNode *gn)
{
  GtFeatureNode *fn = gt_feature_node_cast(gn);
  return fn->seqid;
}

static GtRange feature_node_get_range(GtGenomeNode *gn)
{
  GtFeatureNode *fn = gt_feature_node_cast(gn);
  return fn->range;
}

static void feature_node_set_range(GtGenomeNode *gn, const GtRange *range)
{
  GtFeatureNode *fn = gt_feature_node_cast(gn);
  fn->range = *range;
  if (fn->observer && fn->observer->range_changed)
    fn->observer->range_changed(fn, &(fn->range), fn->observer->data);
}

static void feature_node_change_seqid(GtGenomeNode *gn, GtStr *seqid)
{
  GtFeatureNode *fn = gt_feature_node_cast(gn);
  gt_assert(fn && seqid);
  gt_str_delete(fn->seqid);
  fn->seqid = gt_str_ref(seqid);
}

void gt_feature_node_set_source(GtFeatureNode *fn, GtStr *source)
{
  gt_assert(fn && source);
  if (fn->source)
    gt_str_delete(fn->source);
  fn->source = gt_str_ref(source);
  if (fn->observer && fn->observer->source_changed)
    fn->observer->source_changed(fn, source, fn->observer->data);
}

void gt_feature_node_set_phase(GtFeatureNode *fn, GtPhase phase)
{
  gt_assert(fn);
  fn->bit_field &= ~(PHASE_MASK << PHASE_OFFSET);
  fn->bit_field |= phase << PHASE_OFFSET;
  if (fn->observer && fn->observer->phase_changed)
    fn->observer->phase_changed(fn, phase, fn->observer->data);
}

static int feature_node_accept(GtGenomeNode *gn, GtNodeVisitor *nv,
                               GtError *err)
{
  GtFeatureNode *fn;
  gt_error_check(err);
  fn = gt_feature_node_cast(gn);
  return gt_node_visitor_visit_feature_node(nv, fn, err);
}

void gt_feature_node_set_observer(GtFeatureNode *fn, GtFeatureNodeObserver *o)
{
  gt_assert(fn && o);
  if (fn->observer)
    gt_feature_node_observer_delete(fn->observer);
  fn->observer = gt_feature_node_observer_ref(o);
}

void gt_feature_node_unset_observer(GtFeatureNode *fn)
{
  gt_assert(fn);
  if (fn->observer)
    gt_feature_node_observer_delete(fn->observer);
  fn->observer = NULL;
}

const GtGenomeNodeClass* gt_feature_node_class()
{
  static const GtGenomeNodeClass *gnc = NULL;
  gt_class_alloc_lock_enter();
  if (!gnc) {
    gnc = gt_genome_node_class_new(sizeof (GtFeatureNode),
                                   feature_node_free,
                                   feature_node_get_seqid,
                                   feature_node_get_seqid,
                                   feature_node_get_range,
                                   feature_node_set_range,
                                   feature_node_change_seqid,
                                   feature_node_accept);
  }
  gt_class_alloc_lock_leave();
  return gnc;
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
                                  GtUword start, GtUword end,
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
  fn->score       = GT_UNDEF_FLOAT;
  fn->range.start = start;
  fn->range.end   = end;
  fn->representative = NULL;
  fn->attributes  = NULL;
  fn->bit_field   = 0;
  fn->bit_field  |= strand << STRAND_OFFSET;
  fn->children    = NULL; /* the children list is create on demand */
  fn->observer    = NULL;
  gt_feature_node_set_phase(fn, GT_PHASE_UNDEFINED);
  set_transcriptfeaturetype(fn, TRANSCRIPT_FEATURE_TYPE_UNDETERMINED);
  set_tree_status(&fn->bit_field, IS_TREE);
  /* the DFS status is set to DFS_WHITE already */
  fn->representative = NULL;
  return gn;
}

GtGenomeNode* gt_feature_node_new_pseudo(GtStr *seqid, GtUword start,
                                         GtUword end, GtStrand strand)
{
  GtFeatureNode *pf;
  GtGenomeNode *pn;
  gt_assert(seqid);
  gt_assert(start <= end);
  pn = gt_feature_node_new(seqid, "pseudo", start, end, strand);
  pf = gt_feature_node_cast(pn);
  pf->type = NULL; /* pseudo features do not have a type */
  pf->bit_field |= 1 << PSEUDO_FEATURE_OFFSET;
  return pn;
}

GtGenomeNode* gt_feature_node_new_pseudo_template(GtFeatureNode *fn)
{
  GtFeatureNode *pf;
  GtGenomeNode *pn;
  GtRange range;
  gt_assert(fn);
  range = feature_node_get_range((GtGenomeNode*) fn),
  pn = gt_feature_node_new_pseudo(feature_node_get_seqid((GtGenomeNode*) fn),
                                  range.start, range.end,
                                  gt_feature_node_get_strand(fn));
  pf = gt_feature_node_cast(pn);
  gt_feature_node_set_source(pf, fn->source);
  return pn;
}

GtGenomeNode* gt_feature_node_new_standard_gene(void)
{
  GtGenomeNode *fn, *child, *grand;
  GtStr *seqid;
  seqid = gt_str_new_cstr("ctg123");

  /* gene */
  fn = gt_feature_node_new(seqid, gt_ft_gene, 1000, 9000, GT_STRAND_FORWARD);

  /* TF binding site */
  child = gt_feature_node_new(seqid, gt_ft_TF_binding_site, 1000, 1012,
                                GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) fn, (GtFeatureNode*) child);

  /* first mRNA */
  child = gt_feature_node_new(seqid, gt_ft_mRNA, 1050, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) fn, (GtFeatureNode*) child);
  grand = gt_feature_node_new(seqid, gt_ft_exon, 1050, 1500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grand);
  grand = gt_feature_node_new(seqid, gt_ft_exon, 3000, 3902, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grand);
  grand = gt_feature_node_new(seqid, gt_ft_exon, 5000, 5500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grand);
  grand = gt_feature_node_new(seqid, gt_ft_exon, 7000, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grand);

  /* second mRNA */
  child = gt_feature_node_new(seqid, gt_ft_mRNA, 1050, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) fn, (GtFeatureNode*) child);
  grand = gt_feature_node_new(seqid, gt_ft_exon, 1050, 1500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grand);
  grand = gt_feature_node_new(seqid, gt_ft_exon, 5000, 5500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grand);
  grand = gt_feature_node_new(seqid, gt_ft_exon, 7000, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grand);

  /* third mRNA */
  child = gt_feature_node_new(seqid, gt_ft_mRNA, 1300, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) fn, (GtFeatureNode*) child);
  grand = gt_feature_node_new(seqid, gt_ft_exon, 1300, 1500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grand);
  grand = gt_feature_node_new(seqid, gt_ft_exon, 3000, 3902, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grand);
  grand = gt_feature_node_new(seqid, gt_ft_exon, 5000, 5500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grand);
  grand = gt_feature_node_new(seqid, gt_ft_exon, 7000, 9000,GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) child, (GtFeatureNode*) grand);

  gt_str_delete(seqid);
  return fn;
}

GtFeatureNode* gt_feature_node_clone(const GtFeatureNode *template)
{
  GtFeatureNode *fn;
  GtStr *seqid;
  const char *type;
  GtUword i, start, end;
  GtStrand strand;
  GtStrArray *attributes;
  gt_assert(template && !gt_feature_node_has_children(template));
  seqid = gt_genome_node_get_seqid((GtGenomeNode*) template);
  type = gt_feature_node_get_type(template);
  start = gt_genome_node_get_start((GtGenomeNode*) template);
  end = gt_genome_node_get_end((GtGenomeNode*) template);
  strand = gt_feature_node_get_strand(template);
  fn = (GtFeatureNode*) gt_feature_node_new(seqid, type, start, end, strand);
  gt_genome_node_set_origin((GtGenomeNode*) fn,
                            template->parent_instance.filename,
                            template->parent_instance.line_number);
  if (gt_feature_node_has_source(template))
    gt_feature_node_set_source(fn, template->source);
  if (gt_feature_node_score_is_defined(template))
    gt_feature_node_set_score(fn, template->score);
  attributes = gt_feature_node_get_attribute_list(template);
  for (i = 0; i < gt_str_array_size(attributes); i++) {
    const char *tag, *value;
    tag = gt_str_array_get(attributes, i);
    value = gt_feature_node_get_attribute(template, tag);
    gt_feature_node_set_attribute(fn, tag, value);
  }
  gt_str_array_delete(attributes);
  return fn;
}

const char* gt_feature_node_get_source(const GtFeatureNode *fn)
{
  gt_assert(fn);
  return fn->source ? gt_str_get(fn->source) : ".";
}

bool gt_feature_node_has_source(const GtFeatureNode *fn)
{
  gt_assert(fn);
  if (!fn->source || !strcmp(gt_str_get(fn->source), "."))
    return false;
  return true;
}

const char* gt_feature_node_get_type(const GtFeatureNode *fn)
{
  gt_assert(fn);
  return fn->type;
}

void gt_feature_node_set_type(GtFeatureNode *fn, const char *type)
{
  gt_assert(fn && type);
  fn->type = gt_symbol(type);
  if (fn->observer && fn->observer->type_changed)
    fn->observer->type_changed(fn, fn->type, fn->observer->data);
}

bool gt_feature_node_has_type(GtFeatureNode *fn, const char *type)
{
  gt_assert(fn && type);
  if (!fn->type) /* <fn> is pseudo-node */
    return false;
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
    gt_assert(!((fn->bit_field >> PSEUDO_FEATURE_OFFSET) &
                   PSEUDO_FEATURE_MASK));
    return true;
  }
  return false;
}

bool gt_feature_node_is_pseudo(const GtFeatureNode *fn)
{
  gt_assert(fn);
  if ((fn->bit_field >> PSEUDO_FEATURE_OFFSET) & PSEUDO_FEATURE_MASK) {
    gt_assert(!((fn->bit_field >> MULTI_FEATURE_OFFSET) & MULTI_FEATURE_MASK));
    return true;
  }
  return false;
}

static void feature_node_set_multi(GtFeatureNode *fn)
{
  gt_assert(fn);
  fn->bit_field |= 1 << MULTI_FEATURE_OFFSET;
}

void gt_feature_node_make_multi_representative(GtFeatureNode *fn)
{
  gt_assert(fn);
  feature_node_set_multi(fn);
  if (fn->observer && fn->observer->multi_changed) {
    fn->observer->multi_changed(fn, gt_feature_node_is_multi(fn),
                                fn->representative, fn->observer->data);
  }
}

void gt_feature_node_set_multi_representative(GtFeatureNode *fn,
                                              GtFeatureNode *rep)
{
  gt_assert(fn);
  gt_assert(rep && gt_feature_node_is_multi(rep));
  feature_node_set_multi(fn);
  fn->representative = rep;
  if (fn->observer && fn->observer->multi_changed) {
    fn->observer->multi_changed(fn, gt_feature_node_is_multi(fn),
                                fn->representative, fn->observer->data);
  }
}

void gt_feature_node_unset_multi(GtFeatureNode *fn)
{
  gt_assert(fn);
  fn->bit_field &= ~(1 << MULTI_FEATURE_OFFSET);
  fn->representative = NULL;
}

GtFeatureNode* gt_feature_node_get_multi_representative(GtFeatureNode *fn)
{
  gt_assert(fn && gt_feature_node_is_multi(fn) &&
            !gt_feature_node_is_pseudo(fn));
  if (fn->representative) {
    gt_assert(gt_feature_node_is_multi(fn->representative));
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

void gt_feature_node_get_score_p(const GtFeatureNode *fn, float *val)
{
  gt_assert(fn);
  gt_assert(gt_feature_node_score_is_defined(fn));
  *val = fn->score;
}

GtStrand gt_feature_node_get_strand(const GtFeatureNode *fn)
{
  gt_assert(fn);
  return (fn->bit_field >> STRAND_OFFSET) & STRAND_MASK;
}

void gt_feature_node_set_strand(GtFeatureNode *fn, GtStrand strand)
{
  gt_assert(fn);
  fn->bit_field &= ~(STRAND_MASK << STRAND_OFFSET);
  fn->bit_field |= strand << STRAND_OFFSET;
  if (fn->observer && fn->observer->strand_changed)
    fn->observer->strand_changed(fn, strand, fn->observer->data);
}

GtPhase gt_feature_node_get_phase(const GtFeatureNode *fn)
{
  gt_assert(fn);
  return (fn->bit_field >> PHASE_OFFSET) & PHASE_MASK;
}

static int feature_node_save_exon(GtFeatureNode *fn, void *data,
                                  GT_UNUSED GtError *err)
{
  GtArray *exon_features = (GtArray*) data;
  gt_error_check(err);
  gt_assert(fn && exon_features);
  if (gt_feature_node_has_type(fn, gt_ft_exon)) {
    gt_array_add(exon_features, fn);
  }
  return 0;
}

void gt_feature_node_get_exons(GtFeatureNode *fn, GtArray *exon_features)
{
  GT_UNUSED int had_err;
  gt_assert(fn && exon_features && !gt_array_size(exon_features));
  had_err = gt_feature_node_traverse_children(fn, exon_features,
                                              feature_node_save_exon, false,
                                              NULL);
  gt_assert(!had_err); /* feature_node_save_exon() is sane */
}

static int save_exons_and_cds(GtFeatureNode *fn, void *data,
                              GT_UNUSED GtError *err)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  gt_error_check(err);
  gt_assert(fn && info);
  if (gt_feature_node_has_type(fn, gt_ft_exon))
    gt_array_add(info->exon_features, fn);
  else if (gt_feature_node_has_type(fn, gt_ft_CDS))
    gt_array_add(info->cds_features, fn);
  return 0;
}

static void set_transcript_types(GtArray *features, GtStrand parent_strand)
{
  GtFeatureNode *fn;
  GtUword i;
  gt_assert(features);
  if (gt_array_size(features)) {
    if (gt_array_size(features) == 1) {
      fn = *(GtFeatureNode**) gt_array_get_first(features);
      set_transcriptfeaturetype(fn, TRANSCRIPT_FEATURE_TYPE_SINGLE);
    }
    else {
      fn = *(GtFeatureNode**) gt_array_get_first(features);
      set_transcriptfeaturetype(fn, parent_strand != GT_STRAND_REVERSE
                                    ? TRANSCRIPT_FEATURE_TYPE_INITIAL
                                    : TRANSCRIPT_FEATURE_TYPE_TERMINAL);
      for (i = 1; i < gt_array_size(features) - 1; i++) {
        fn = *(GtFeatureNode**) gt_array_get(features, i);
        set_transcriptfeaturetype(fn, TRANSCRIPT_FEATURE_TYPE_INTERNAL);
      }
      fn = *(GtFeatureNode**) gt_array_get_last(features);
      set_transcriptfeaturetype(fn, parent_strand != GT_STRAND_REVERSE
                                    ? TRANSCRIPT_FEATURE_TYPE_TERMINAL
                                    : TRANSCRIPT_FEATURE_TYPE_INITIAL);
    }
  }
}

static int determine_transcripttypes(GtFeatureNode *fn, void *data,
                                     GT_UNUSED GtError *err)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  GT_UNUSED int had_err;
  gt_error_check(err);
  gt_assert(fn && info);
  /* reset exon_features and cds_features */
  gt_array_reset(info->exon_features);
  gt_array_reset(info->cds_features);
  /* collect all direct children exons */
  had_err = gt_feature_node_traverse_direct_children(fn, info,
                                                     save_exons_and_cds, NULL);
  gt_assert(!had_err); /* cannot happen, because save_exons_and_cds() is sane */
  /* set transcript feature type, if necessary */
  set_transcript_types(info->exon_features, gt_feature_node_get_strand(fn));
  set_transcript_types(info->cds_features, gt_feature_node_get_strand(fn));
  return 0;
}

void gt_feature_node_determine_transcripttypes(GtFeatureNode *fn)
{
  SaveExonAndCDSInfo info;
  GT_UNUSED int had_err;
  gt_assert(fn);
  info.exon_features = gt_array_new(sizeof (GtFeatureNode*));
  info.cds_features = gt_array_new(sizeof (GtFeatureNode*));
  had_err = gt_feature_node_traverse_children(fn, &info,
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

void gt_feature_node_set_end(GtFeatureNode *fn, GtUword end)
{
  gt_assert(fn && fn->range.start <= end);
  fn->range.end = end;
  if (fn->observer && fn->observer->range_changed)
    fn->observer->range_changed(fn, &(fn->range), fn->observer->data);
}

void gt_feature_node_set_score(GtFeatureNode *fn, float score)
{
  gt_assert(fn);
  fn->bit_field |= 1 << SCORE_IS_DEFINED_OFFSET;
  fn->score = score;
  if (fn->observer && fn->observer->score_changed)
    fn->observer->score_changed(fn, score, fn->observer->data);
}

void gt_feature_node_set_score_p(GtFeatureNode *fn, float* score)
{
  gt_assert(fn);
  fn->bit_field |= 1 << SCORE_IS_DEFINED_OFFSET;
  fn->score = *score;
  if (fn->observer && fn->observer->score_changed)
    fn->observer->score_changed(fn, *score, fn->observer->data);
}

void gt_feature_node_unset_score(GtFeatureNode *fn)
{
  gt_assert(fn);
  fn->bit_field &= ~(1 << SCORE_IS_DEFINED_OFFSET);
  fn->score = GT_UNDEF_FLOAT;
  if (fn->observer && fn->observer->score_changed)
    fn->observer->score_changed(fn, fn->score, fn->observer->data);
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
  if (fn->observer && fn->observer->attribute_changed) {
    fn->observer->attribute_changed(fn, true, attr_name, attr_value,
                                    fn->observer->data);
  }
}

void gt_feature_node_set_attribute(GtFeatureNode *fn,
                                   const char *attr_name,
                                   const char *attr_value)
{
  gt_assert(fn && attr_name && attr_value);
  gt_assert(strlen(attr_name)); /* attribute name cannot be empty */
  gt_assert(strlen(attr_value)); /* attribute value cannot be empty */
  if (!fn->attributes)
    fn->attributes = gt_tag_value_map_new(attr_name, attr_value);
  else
    gt_tag_value_map_set(&fn->attributes, attr_name, attr_value);
  if (fn->observer && fn->observer->attribute_changed) {
    fn->observer->attribute_changed(fn, false, attr_name, attr_value,
                                    fn->observer->data);
  }
}

void gt_feature_node_remove_attribute(GtFeatureNode *fn,
                                      const char *attr_name)
{
  gt_assert(fn && attr_name);
  gt_assert(strlen(attr_name)); /* attribute name cannot be empty */
  gt_assert(fn->attributes); /* attribute list must exist already */
  gt_tag_value_map_remove(&fn->attributes, attr_name);
  if (fn->observer && fn->observer->attribute_deleted) {
    fn->observer->attribute_deleted(fn, attr_name, fn->observer->data);
  }
}

void gt_feature_node_foreach_attribute(GtFeatureNode *fn,
                                       GtFeatureNodeAttributeIterFunc iterfunc,
                                       void *data)
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
    GtUword i = 0;
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
  static const char *gfts[] = { gt_ft_CDS, NULL };
  return feature_node_has_gft(fn, gfts);
}

bool gt_feature_node_has_splice_site(const GtFeatureNode *fn)
{
  static const char *gfts[] = { gt_ft_five_prime_cis_splice_site,
                                gt_ft_five_prime_splice_site,
                                gt_ft_three_prime_cis_splice_site,
                                gt_ft_three_prime_splice_site,
                                NULL };
  return feature_node_has_gft(fn, gfts);
}

double gt_feature_node_average_splice_site_prob(const GtFeatureNode *fn,
                                                GtUword *num_ss)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *child;
  GtUword num_of_splice_sites = 0;
  double averagessp = 0.0;
  gt_assert(fn);
  fni = gt_feature_node_iterator_new(fn);
  while ((child = gt_feature_node_iterator_next(fni))) {
    if (gt_feature_node_has_type(child, gt_ft_five_prime_cis_splice_site) ||
        gt_feature_node_has_type(child, gt_ft_five_prime_splice_site) ||
        gt_feature_node_has_type(child, gt_ft_three_prime_cis_splice_site) ||
        gt_feature_node_has_type(child, gt_ft_three_prime_splice_site)) {
      averagessp += gt_feature_node_get_score((GtFeatureNode*) child);
      num_of_splice_sites++;
    }
  }
  gt_feature_node_iterator_delete(fni);
  if (num_of_splice_sites)
    averagessp /= num_of_splice_sites;
  if (num_ss)
    *num_ss = num_of_splice_sites;
  return averagessp;
}

bool gt_feature_node_is_similar(const GtFeatureNode *fn_a,
                                const GtFeatureNode *fn_b)
{
  GtRange range_a, range_b;
  gt_assert(fn_a && fn_b);
  range_a = feature_node_get_range((GtGenomeNode*) fn_a);
  range_b = feature_node_get_range((GtGenomeNode*) fn_b);
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
  fn = gt_feature_node_new(seqid, gt_ft_gene, 1, 1000, GT_STRAND_FORWARD);

  gt_ensure(!gt_feature_node_score_is_defined((GtFeatureNode*) fn));

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

int gt_feature_node_traverse_children(GtFeatureNode *feature_node, void *data,
                                      GtFeatureNodeTraverseFunc traverse,
                                      bool traverse_only_once, GtError *err)
{
  GtArray *node_stack = NULL, *list_of_children;
  GtFeatureNode *fn, *fn_ref, *child_feature;
  GtDlistelem *dlistelem;
  GtUword i;
  GtHashtable *traversed_nodes = NULL;
  bool has_node_with_multiple_parents = false;
  int had_err = 0;

  if (!feature_node)
    return 0;

  /* create additional reference to <feature_node> (necessary if feature_node is
     freed by <traverse>) */
  fn_ref = (GtFeatureNode*) gt_genome_node_ref((GtGenomeNode*) feature_node);

  node_stack = gt_array_new(sizeof (GtFeatureNode*));
  if (gt_feature_node_is_pseudo(feature_node)) {
    /* add the children backwards to traverse in order */
    for (dlistelem = gt_dlist_last(feature_node->children);
         dlistelem != NULL;
         dlistelem = gt_dlistelem_previous(dlistelem)) {
      child_feature = (GtFeatureNode*) gt_dlistelem_get_data(dlistelem);
      gt_array_add(node_stack, child_feature);
    }
  }
  else
    gt_array_add(node_stack, feature_node);
  gt_assert(gt_array_size(node_stack));
  list_of_children = gt_array_new(sizeof (GtFeatureNode*));

  if (traverse_only_once) {
    static const HashElemInfo node_hashtype
      = { gt_ht_ptr_elem_hash, { NULL }, sizeof (GtFeatureNode *),
          gt_ht_ptr_elem_cmp, NULL, NULL };
    traversed_nodes = gt_hashtable_new(node_hashtype);
  }

  while (gt_array_size(node_stack)) {
    fn = *(GtFeatureNode**) gt_array_pop(node_stack);
    gt_array_reset(list_of_children);
    if (fn->children) {
      /* a backup of the children array is necessary if traverse() frees the
         node */
      for (dlistelem = gt_dlist_first(fn->children); dlistelem != NULL;
           dlistelem = gt_dlistelem_next(dlistelem)) {
        child_feature = (GtFeatureNode*) gt_dlistelem_get_data(dlistelem);
        gt_array_add(list_of_children, child_feature);
      }
    }
    /* store the implications of <fn> to the tree status of <feature_node> */
    if (multiple_parents(fn->bit_field))
      has_node_with_multiple_parents = true;
    /* call traverse function */
    if (traverse) {
      if ((had_err = traverse(fn, data, err)))
        break;
    }
    for (i = 0; i < gt_array_size(list_of_children); i++) {
      /* we go backwards to traverse in order */
      child_feature = *(GtFeatureNode**) gt_array_get(list_of_children,
                                       gt_array_size(list_of_children) - i - 1);
      if (!traverse_only_once ||
          !gt_hashtable_get(traversed_nodes, &child_feature)) {
        /* feature has not been traversed or has to be traversed multiple
           times */
        gt_array_add(node_stack, child_feature);
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
  gt_genome_node_delete((GtGenomeNode*) fn_ref);
  if (traverse_only_once)
    gt_hashtable_delete(traversed_nodes);
  gt_array_delete(list_of_children);
  gt_array_delete(node_stack);

  return had_err;
}

DFSStatus feature_node_get_dfs_status(const GtFeatureNode *fn)
{
  gt_assert(fn);
  return (fn->bit_field >> DFS_STATUS_OFFSET) & DFS_STATUS_MASK;
}

void feature_node_set_dfs_status(GtFeatureNode *fn, DFSStatus status)
{
  gt_assert(fn);
  fn->bit_field &= ~(DFS_STATUS_MASK << DFS_STATUS_OFFSET);
  fn->bit_field |= status << DFS_STATUS_OFFSET;
}

static void dfs_visit(GtFeatureNode *u, GtArray *toplist)
{
  GtDlistelem *dlistelem;
  gt_assert(u && toplist);
  feature_node_set_dfs_status(u, DFS_GRAY);
  if (u->children) {
    for (dlistelem = gt_dlist_last(u->children);
         dlistelem != NULL;
         dlistelem = gt_dlistelem_previous(dlistelem)) {
      GtFeatureNode *v = (GtFeatureNode*) gt_dlistelem_get_data(dlistelem);
      if (feature_node_get_dfs_status(v) == DFS_WHITE)
        dfs_visit(v, toplist);
    }
  }
  feature_node_set_dfs_status(u, DFS_BLACK);
  if (!gt_feature_node_is_pseudo(u))
    gt_array_add(toplist, u);
}

/* Implements topologically sorted depth-first search on directed graphs.
   For a description see for example chapter 23 of the book:

   T.H. Cormen, C.E. Leiserson and R.L. Rivest. __Introduction to Algorithms__.
   MIT Press: Cambridge, MA, 1990. */
int gt_feature_node_traverse_children_top(GtFeatureNode *feature_node,
                                          void *data,
                                          GtFeatureNodeTraverseFunc traverse,
                                          GtError *err)
{
  GtArray *toplist;
  int had_err = 0;

  if (!feature_node)
    return 0;

  /* this is the first traversal */
  gt_assert(feature_node_get_dfs_status(feature_node) == DFS_WHITE);

  /* initialization */
  toplist = gt_array_new(sizeof (GtFeatureNode*));

  /* process queue */
  dfs_visit(feature_node, toplist);

  /* process topologically sorted feature node list */
  if (traverse) {
    while (gt_array_size(toplist)) {
      GtFeatureNode *fn = *(GtFeatureNode**) gt_array_pop(toplist);
      if ((had_err = traverse(fn, data, err)))
        break;
    }
  }

  /* free */
  gt_array_delete(toplist);

  return had_err;
}

static int count_types(GtFeatureNode *fn, void *data, GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtTypeTraverseInfo *traverseinfo = (GtTypeTraverseInfo*) data;
  gt_error_check(err);
  if (strcmp(traverseinfo->type, gt_feature_node_get_type(fn)) == 0)
    traverseinfo->number++;
  return had_err;
}

int gt_feature_node_traverse_direct_children(GtFeatureNode *fn,
                                             void *traverse_func_data,
                                             GtFeatureNodeTraverseFunc traverse,
                                             GtError *err)
{
  GtDlistelem *dlistelem;
  int had_err = 0;
  gt_error_check(err);
  if (!fn || !traverse)
    return 0;
  if (fn->children) {
    for (dlistelem = gt_dlist_first(fn->children); dlistelem != NULL;
         dlistelem = gt_dlistelem_next(dlistelem)) {
      had_err = traverse((GtFeatureNode*) gt_dlistelem_get_data(dlistelem),
                          traverse_func_data, err);
      if (had_err)
        break;
    }
  }
  return had_err;
}

GtUword gt_feature_node_number_of_children(const GtFeatureNode *fn)
{
  gt_assert(fn);
  return gt_dlist_size(fn->children);
}

GtUword gt_feature_node_number_of_children_of_type(const GtFeatureNode
                                                         *parent,
                                                         const GtFeatureNode
                                                         *node)
{
  GT_UNUSED int had_err = 0;
  GtTypeTraverseInfo traverseinfo;
  gt_assert(parent && node);
  traverseinfo.type = gt_feature_node_get_type(node);
  traverseinfo.number = 0;
  had_err = gt_feature_node_traverse_direct_children((GtFeatureNode*) parent,
                                                     &traverseinfo, count_types,
                                                     NULL);
  return traverseinfo.number;
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
  gt_dlist_add(parent->children, child); /* XXX: check for cycles */
  /* update tree status of <parent> */
  set_tree_status(&parent->bit_field, TREE_STATUS_UNDETERMINED);
  /* update parent info of <child> */
  add_parent(&child->bit_field);
  if (parent->observer && parent->observer->child_added)
    parent->observer->child_added(parent, child, parent->observer->data);
}

typedef struct {
  GtFeatureNode *leaf,
                *parent;
  bool deleted;
} GtFeatureNodeLeafDeleteInfo;

static int remove_leaf(GtFeatureNode *node, void *data, GT_UNUSED GtError *err)
{
  GtFeatureNode *child;
  GtFeatureNodeLeafDeleteInfo* info = (GtFeatureNodeLeafDeleteInfo*) data;
  GtDlistelem *dlistelem;
  gt_error_check(err);
  if (node != info->leaf && node->children) {
    for (dlistelem = gt_dlist_first(node->children); dlistelem != NULL;
         dlistelem = gt_dlistelem_next(dlistelem)) {
      child = (GtFeatureNode*) gt_dlistelem_get_data(dlistelem);
      if (child == info->leaf) {
        gt_dlist_remove(node->children, dlistelem);
        gt_genome_node_delete((GtGenomeNode*) child);
        info->deleted = true;
        break;
      }
    }
  }
  return 0;
}

void gt_feature_node_remove_leaf(GtFeatureNode *tree, GtFeatureNode *leafn)
{
  GT_UNUSED int had_err;
  GtFeatureNodeLeafDeleteInfo info;
  gt_assert(tree && leafn);

  info.leaf = leafn;
  info.parent = tree;
  info.deleted = false;

  /* ref child node to enable traversal */
  gt_genome_node_ref((GtGenomeNode*) leafn);
  gt_assert(gt_feature_node_number_of_children(leafn) == 0);
  had_err = gt_feature_node_traverse_children(tree, &info, remove_leaf, true,
                                              NULL);
  /* unref child node, traversal done */
  gt_genome_node_delete((GtGenomeNode*) leafn);

  if (info.deleted)
    set_tree_status(&tree->bit_field, TREE_STATUS_UNDETERMINED);

  gt_assert(!had_err); /* cannot happen, remove_leaf() is sane */
}

void gt_feature_node_mark(GtFeatureNode *fn)
{
  gt_assert(fn);
  fn->bit_field |= 1;
  if (fn->observer && fn->observer->mark_changed)
    fn->observer->mark_changed(fn, true, fn->observer->data);
}

void gt_feature_node_unmark(GtFeatureNode *fn)
{
  gt_assert(fn);
  if (gt_feature_node_is_marked(fn))
    fn->bit_field &= ~(1);
  if (fn->observer && fn->observer->mark_changed)
    fn->observer->mark_changed(fn, false, fn->observer->data);
}

bool gt_feature_node_is_marked(const GtFeatureNode *fn)
{
  gt_assert(fn);
  return fn->bit_field & 1 ? true : false;
}

static int check_marked_status(GtFeatureNode *fn, void *data,
                               GT_UNUSED GtError *err)
{
  bool *marked = data;
  if (gt_feature_node_is_marked(fn))
    *marked = true;
  return 0;
}

bool gt_feature_node_contains_marked(GtFeatureNode *fn)
{
  bool contains_marked = false;
  GT_UNUSED int rval;
  gt_assert(fn);
  rval = gt_feature_node_traverse_children(fn, &contains_marked,
                                           check_marked_status, true, NULL);
  gt_assert(!rval); /* check_marked_status() is sane */
  return contains_marked;
}

bool gt_feature_node_has_children(const GtFeatureNode *fn)
{
  gt_assert(fn);
  if (!fn->children || gt_dlist_size(fn->children) == 0)
    return false;
  return true;
}

bool gt_feature_node_direct_children_do_not_overlap_generic(GtFeatureNode
                                                            *parent,
                                                            GtFeatureNode
                                                            *child)
{
  GtArray *children_ranges;
  GtDlistelem *dlistelem;
  GtFeatureNode *child_fn;
  GtRange range;
  bool rval;

  gt_assert(parent);

  if (!parent->children)
    return true;

  /* get children ranges */
  children_ranges = gt_array_new(sizeof (GtRange));
  gt_assert(parent->children);
  for (dlistelem = gt_dlist_first(parent->children); dlistelem != NULL;
       dlistelem = gt_dlistelem_next(dlistelem)) {
    if (!child ||
        ((child_fn =
            gt_feature_node_try_cast(gt_dlistelem_get_data(dlistelem))) &&
         gt_feature_node_get_type(child) ==
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

bool gt_feature_node_direct_children_do_not_overlap(GtFeatureNode *fn)
{
  return gt_feature_node_direct_children_do_not_overlap_generic(fn, NULL);
}

bool gt_feature_node_direct_children_do_not_overlap_st(GtFeatureNode *parent,
                                                       GtFeatureNode *child)
{
  return gt_feature_node_direct_children_do_not_overlap_generic(parent, child);
}

bool gt_feature_node_is_tree(GtFeatureNode *fn)
{
  bool status = false;
  gt_assert(fn);
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

bool gt_feature_node_overlaps_nodes(GtFeatureNode *fn, GtArray *nodes)
{
  return gt_feature_node_overlaps_nodes_mark(fn, nodes, NULL);
}

bool gt_feature_node_overlaps_nodes_mark(GtFeatureNode *fn, GtArray *nodes,
                                         GtBittab *b)
{
  GtUword i;
  GtGenomeNode *node;
  GtRange fn_range, node_range;
  bool rval = false;
#ifndef NDEBUG
  GtStr *fn_id;
  gt_assert(fn && nodes);
  gt_assert(!b || gt_bittab_size(b) == gt_array_size(nodes));
  fn_id = gt_genome_node_get_idstr((GtGenomeNode*) fn);
#endif
  fn_range = gt_genome_node_get_range((GtGenomeNode*) fn);

  for (i = 0; i < gt_array_size(nodes); i++) {
    node = *(GtGenomeNode**) gt_array_get(nodes, i);
    gt_assert(!gt_str_cmp(fn_id, gt_genome_node_get_idstr(node)));
    node_range = gt_genome_node_get_range(node);
    if (gt_range_overlap(&fn_range, &node_range)) {
      rval = true;
      if (b)
        gt_bittab_set_bit(b, i);
      else
        break;
    }
  }
  return rval;
}

GtFeatureNode* gt_feature_node_try_cast(GtGenomeNode *gn)
{
  return gt_genome_node_try_cast(gt_feature_node_class(), gn);
}

GtFeatureNode* gt_feature_node_cast(GtGenomeNode *gn)
{
  return gt_genome_node_cast(gt_feature_node_class(), gn);
}
