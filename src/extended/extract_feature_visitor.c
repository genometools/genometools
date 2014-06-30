/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
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

#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/codon_iterator_simple_api.h"
#include "core/fasta.h"
#include "core/symbol_api.h"
#include "core/translator.h"
#include "extended/extract_feature_sequence.h"
#include "extended/extract_feature_visitor.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"

struct GtExtractFeatureVisitor {
  const GtNodeVisitor parent_instance;
  const char *type;
  bool join,
       translate,
       seqid,
       target,
       retain_ids;
  GtUword fastaseq_counter,
                width;
  GtRegionMapping *region_mapping;
  GtFile *outfp;
};

#define gt_extract_feature_visitor_cast(GV)\
        gt_node_visitor_cast(gt_extract_feature_visitor_class(), GV)

static void extract_feature_visitor_free(GtNodeVisitor *nv)
{
  GtExtractFeatureVisitor *efv = gt_extract_feature_visitor_cast(nv);
  gt_assert(efv);
  gt_region_mapping_delete(efv->region_mapping);
}

static void construct_description(GtStr *description, const char *type,
                                  GtUword counter, bool join,
                                  bool translate, GtStr *seqid,
                                  GtStrArray *target_ids)
{
  if (!gt_str_length(description)) {
    gt_str_append_cstr(description, type);
    gt_str_append_char(description, '_');
    gt_str_append_ulong(description, counter);
  }
  if (join)
    gt_str_append_cstr(description, " (joined)");
  if (translate)
    gt_str_append_cstr(description, " (translated)");
  if (seqid) {
    gt_assert(gt_str_length(seqid));
    gt_str_append_cstr(description, " [seqid ");
    gt_str_append_str(description, seqid);
    gt_str_append_char(description, ']');
  }
  if (target_ids && gt_str_array_size(target_ids)) {
    GtUword i;
    gt_str_append_cstr(description, " [target IDs ");
    gt_str_append_cstr(description, gt_str_array_get(target_ids, 0));
    for (i = 1; i < gt_str_array_size(target_ids); i++) {
      gt_str_append_char(description, ',');
      gt_str_append_cstr(description, gt_str_array_get(target_ids, i));
    }
    gt_str_append_char(description, ']');
  }
}

static int extract_feature_visitor_feature_node(GtNodeVisitor *nv,
                                                GtFeatureNode *fn, GtError *err)
{
  GtExtractFeatureVisitor *efv;
  GtFeatureNodeIterator *fni;
  GtFeatureNode *child;
  GtStrArray *target_ids = NULL;
  GtStr *seqid = NULL,
        *description,
        *sequence;
  int had_err = 0;
  gt_error_check(err);
  efv = gt_extract_feature_visitor_cast(nv);
  gt_assert(efv->region_mapping);
  fni = gt_feature_node_iterator_new(fn);
  if (efv->target)
    target_ids = gt_str_array_new();
  if (efv->seqid)
    seqid = gt_str_new();
  description = gt_str_new();
  sequence = gt_str_new();
  while (!had_err && (child = gt_feature_node_iterator_next(fni))) {
    if (seqid)
      gt_str_reset(seqid);
    if (target_ids)
      gt_str_array_reset(target_ids);
    if (efv->translate) {
      if (gt_extract_and_translate_feature_sequence(child,
                                                    efv->type,
                                                    efv->region_mapping, NULL,
                                                    sequence, NULL, NULL,
                                                    err)) {
        had_err = -1;
      }
    } else {
      if (gt_extract_feature_sequence(sequence, (GtGenomeNode*) child, efv->type,
                                      efv->join, seqid, target_ids,
                                      efv->region_mapping, err)) {
        had_err = -1;
      }
    }

    if (!had_err && gt_str_length(sequence)) {
      efv->fastaseq_counter++;
      if (efv->retain_ids && gt_feature_node_get_attribute(child, "ID")) {
        gt_assert(!gt_str_length(description));
        gt_str_append_cstr(description, gt_feature_node_get_attribute(child,
                                                                      "ID"));
      }
      construct_description(description, efv->type, efv->fastaseq_counter,
                            efv->join, efv->translate, seqid, target_ids);
      gt_fasta_show_entry(gt_str_get(description), gt_str_get(sequence),
                          gt_str_length(sequence),  efv->width,  efv->outfp);
      gt_str_reset(description);
      gt_str_reset(sequence);
    }

  }
  gt_str_delete(sequence);
  gt_str_delete(description);
  gt_str_delete(seqid);
  gt_str_array_delete(target_ids);
  gt_feature_node_iterator_delete(fni);
  return had_err;
}

const GtNodeVisitorClass* gt_extract_feature_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtExtractFeatureVisitor),
                                    extract_feature_visitor_free,
                                    NULL,
                                    extract_feature_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_extract_feature_visitor_new(GtRegionMapping *rm,
                                              const char *type, bool join,
                                              bool translate, bool seqid,
                                              bool target, GtUword width,
                                              GtFile *outfp)
{
  GtNodeVisitor *nv;
  GtExtractFeatureVisitor *efv;
  gt_assert(rm);
  nv = gt_node_visitor_create(gt_extract_feature_visitor_class());
  efv= gt_extract_feature_visitor_cast(nv);
  efv->type = gt_symbol(type);
  efv->join = join;
  efv->translate = translate;
  efv->seqid = seqid;
  efv->target = target;
  efv->fastaseq_counter = 0;
  efv->region_mapping = rm;
  efv->width = width;
  efv->outfp = outfp;
  /* XXX */
  efv->retain_ids = getenv("GT_RETAINIDS") ? true : false;
  return nv;
}

void gt_extract_feature_visitor_retain_id_attributes(GtExtractFeatureVisitor
                                                                          *efv)
{
  gt_assert(efv);
  efv->retain_ids = true;
}
