/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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
#include "extended/script_wrapper_visitor_api.h"
#include "extended/node_visitor_api.h"

struct GtScriptWrapperVisitor {
  const GtNodeVisitor parent_instance;
  GtScriptWrapperVisitorCommentNodeFunc comment_node_func;
  GtScriptWrapperVisitorFeatureNodeFunc feature_node_func;
  GtScriptWrapperVisitorRegionNodeFunc region_node_func;
  GtScriptWrapperVisitorSequenceNodeFunc sequence_node_func;
  GtScriptWrapperVisitorMetaNodeFunc meta_node_func;
  GtScriptWrapperVisitorEOFNodeFunc eof_node_func;
  GtScriptWrapperVisitorFreeFunc free_func;
};

const GtNodeVisitorClass* gt_script_wrapper_visitor_class(void);

#define gt_script_wrapper_visitor_cast(GV)\
        gt_node_visitor_cast(gt_script_wrapper_visitor_class(), GV)

static void script_wrapper_visitor_free(GtNodeVisitor *nv)
{
  GtScriptWrapperVisitor *swv;
  gt_assert(nv);
  swv  = gt_script_wrapper_visitor_cast(nv);
  if (swv->free_func)
    swv->free_func(NULL);
}

static int script_wrapper_visitor_feature_node(GtNodeVisitor *nv,
                                               GtFeatureNode *fn,
                                               GtError *err)
{
  GtScriptWrapperVisitor *swv;
  int had_err = 0;
  gt_error_check(err);
  swv = gt_script_wrapper_visitor_cast(nv);
  if (swv->feature_node_func)
    had_err = swv->feature_node_func(fn, err);
  return had_err;
}

static int script_wrapper_visitor_region_node(GtNodeVisitor *nv,
                                              GtRegionNode *rn,
                                              GtError *err)
{
  GtScriptWrapperVisitor *swv;
  int had_err = 0;
  gt_error_check(err);
  swv = gt_script_wrapper_visitor_cast(nv);
  if (swv->region_node_func)
    had_err = swv->region_node_func(rn, err);
  return had_err;
}

static int script_wrapper_visitor_comment_node(GtNodeVisitor *nv,
                                               GtCommentNode *cn,
                                               GtError *err)
{
  GtScriptWrapperVisitor *swv;
  int had_err = 0;
  gt_error_check(err);
  swv = gt_script_wrapper_visitor_cast(nv);
  if (swv->comment_node_func)
    had_err = swv->comment_node_func(cn, err);
  return had_err;
}

static int script_wrapper_visitor_sequence_node(GtNodeVisitor *nv,
                                                GtSequenceNode *sn,
                                                GtError *err)
{
  GtScriptWrapperVisitor *swv;
  int had_err = 0;
  gt_error_check(err);
  swv = gt_script_wrapper_visitor_cast(nv);
  if (swv->sequence_node_func)
    had_err = swv->sequence_node_func(sn, err);
  return had_err;
}

static int script_wrapper_visitor_meta_node(GtNodeVisitor *nv,
                                            GtMetaNode *rn,
                                            GtError *err)
{
  GtScriptWrapperVisitor *swv;
  int had_err = 0;
  gt_error_check(err);
  swv = gt_script_wrapper_visitor_cast(nv);
  if (swv->meta_node_func)
    had_err = swv->meta_node_func(rn, err);
  return had_err;
}

static int script_wrapper_visitor_eof_node(GtNodeVisitor *nv,
                                           GtEOFNode *rn,
                                           GtError *err)
{
  GtScriptWrapperVisitor *swv;
  int had_err = 0;
  gt_error_check(err);
  swv = gt_script_wrapper_visitor_cast(nv);
  if (swv->eof_node_func)
    had_err = swv->eof_node_func(rn, err);
  return had_err;
}

const GtNodeVisitorClass* gt_script_wrapper_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtScriptWrapperVisitor),
                                    script_wrapper_visitor_free,
                                    script_wrapper_visitor_comment_node,
                                    script_wrapper_visitor_feature_node,
                                    script_wrapper_visitor_region_node,
                                    script_wrapper_visitor_sequence_node,
                                    script_wrapper_visitor_eof_node);
    gt_node_visitor_class_set_meta_node_func(nvc,
                                             script_wrapper_visitor_meta_node);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor*
gt_script_wrapper_visitor_new(GtScriptWrapperVisitorCommentNodeFunc cn,
                              GtScriptWrapperVisitorFeatureNodeFunc fn,
                              GtScriptWrapperVisitorRegionNodeFunc rn,
                              GtScriptWrapperVisitorSequenceNodeFunc sn,
                              GtScriptWrapperVisitorMetaNodeFunc mn,
                              GtScriptWrapperVisitorEOFNodeFunc en,
                              GtScriptWrapperVisitorFreeFunc free_func)
{
  GtNodeVisitor *nv;
  GtScriptWrapperVisitor *swv;
  nv = gt_node_visitor_create(gt_script_wrapper_visitor_class());

  swv = gt_script_wrapper_visitor_cast(nv);
  swv->comment_node_func = cn;
  swv->feature_node_func = fn;
  swv->region_node_func = rn;
  swv->sequence_node_func = sn;
  swv->meta_node_func = mn;
  swv->eof_node_func = en;
  swv->free_func = free_func;
  return nv;
}
