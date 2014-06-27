/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd

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

#include "core/class_alloc_lock.h"
#include "core/fasta_api.h"
#include "core/queue_api.h"
#include "core/str_api.h"
#include "core/unused_api.h"
#include "extended/sequence_node_out_visitor.h"
#include "extended/node_visitor_api.h"

struct GtSequenceNodeOutVisitor {
  const GtNodeVisitor parent_instance;
  GtQueue *node_buffer;
  GtFile *outfile;
  GtUword width;
  bool keep_sequence_nodes;
};

const GtNodeVisitorClass* gt_sequence_node_out_visitor_class(void);

#define sequence_node_out_visitor_cast(NS)\
        gt_node_visitor_cast(gt_sequence_node_out_visitor_class(), NS)

static void sequence_node_out_visitor_free(GtNodeVisitor *nv)
{
  GtSequenceNodeOutVisitor *v = sequence_node_out_visitor_cast(nv);
  gt_queue_delete(v->node_buffer);
}

static int sequence_node_out_visitor_sequence_node(GtNodeVisitor *nv,
                                                   GtSequenceNode *sn,
                                                   GT_UNUSED GtError *err)
{
  GtSequenceNodeOutVisitor *v = sequence_node_out_visitor_cast(nv);
  GT_UNUSED GtFeatureNode *node;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(v && sn);

  gt_fasta_show_entry(gt_sequence_node_get_description(sn),
                      gt_sequence_node_get_sequence(sn),
                      gt_sequence_node_get_sequence_length(sn),
                      v->width, v->outfile);

  if (v->keep_sequence_nodes)
    gt_queue_add(v->node_buffer, sn);
  else
    gt_genome_node_delete((GtGenomeNode*) sn);

  return had_err;
}

static int sequence_node_out_visitor_comment_node(GtNodeVisitor *nv,
                                                   GtCommentNode *cn,
                                                   GT_UNUSED GtError *err)
{
  GtSequenceNodeOutVisitor *v = sequence_node_out_visitor_cast(nv);
  int had_err = 0;
  gt_error_check(err);
  gt_assert(v && cn);

  gt_queue_add(v->node_buffer, cn);

  return had_err;
}

static int sequence_node_out_visitor_region_node(GtNodeVisitor *nv,
                                                 GtRegionNode *rn,
                                                 GT_UNUSED GtError *err)
{
  GtSequenceNodeOutVisitor *v = sequence_node_out_visitor_cast(nv);
  int had_err = 0;
  gt_error_check(err);
  gt_assert(v && rn);

  gt_queue_add(v->node_buffer, rn);

  return had_err;
}

static int sequence_node_out_visitor_eof_node(GtNodeVisitor *nv,
                                              GtEOFNode *en,
                                              GT_UNUSED GtError *err)
{
  GtSequenceNodeOutVisitor *v = sequence_node_out_visitor_cast(nv);
  int had_err = 0;
  gt_error_check(err);
  gt_assert(v && en);

  gt_queue_add(v->node_buffer, en);

  return had_err;
}

static int sequence_node_out_visitor_feature_node(GtNodeVisitor *nv,
                                                  GtFeatureNode *fn,
                                                  GT_UNUSED GtError *err)
{
  GtSequenceNodeOutVisitor *v = sequence_node_out_visitor_cast(nv);
  int had_err = 0;
  gt_error_check(err);
  gt_assert(v && fn);

  gt_queue_add(v->node_buffer, fn);

  return had_err;
}

static int sequence_node_out_visitor_meta_node(GtNodeVisitor *nv,
                                               GtMetaNode *mn,
                                               GT_UNUSED GtError *err)
{
  GtSequenceNodeOutVisitor *v = sequence_node_out_visitor_cast(nv);
  int had_err = 0;
  gt_error_check(err);
  gt_assert(v && mn);

  gt_queue_add(v->node_buffer, mn);

  return had_err;
}

const GtNodeVisitorClass* gt_sequence_node_out_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtSequenceNodeOutVisitor),
                                    sequence_node_out_visitor_free,
                                    sequence_node_out_visitor_comment_node,
                                    sequence_node_out_visitor_feature_node,
                                    sequence_node_out_visitor_region_node,
                                    sequence_node_out_visitor_sequence_node,
                                    sequence_node_out_visitor_eof_node);
    gt_node_visitor_class_set_meta_node_func(nvc,
                                           sequence_node_out_visitor_meta_node);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_sequence_node_out_visitor_new(GtFile *outfile)
{
  GtNodeVisitor *nv = gt_node_visitor_create(gt_sequence_node_out_visitor_class());
  GtSequenceNodeOutVisitor *v = sequence_node_out_visitor_cast(nv);
  v->outfile = outfile;
  v->width = 80;
  v->keep_sequence_nodes = false;
  v->node_buffer = gt_queue_new();
  return nv;
}

void gt_sequence_node_out_visitor_keep_sequence_nodes(GtSequenceNodeOutVisitor
                                                                             *v)
{
  gt_assert(v);
  v->keep_sequence_nodes = true;
}

GtUword gt_sequence_node_out_visitor_node_buffer_size(GtNodeVisitor *nv)
{
  GtSequenceNodeOutVisitor *v = sequence_node_out_visitor_cast(nv);
  return gt_queue_size(v->node_buffer);
}

GtGenomeNode* gt_sequence_node_out_visitor_get_node(GtNodeVisitor *nv)
{
  GtSequenceNodeOutVisitor *v = sequence_node_out_visitor_cast(nv);
  return gt_queue_get(v->node_buffer);
}
