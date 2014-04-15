/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
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
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/comment_node_api.h"
#include "extended/genome_node_rep.h"

struct GtCommentNode
{
  const GtGenomeNode parent_instance;
  char *comment;
  GtStr *comment_str; /* used in gt_comment_node_get_idstr() */
};

static void comment_node_free(GtGenomeNode *gn)
{
  GtCommentNode *c = gt_comment_node_cast(gn);
  gt_assert(c && c->comment);
  gt_free(c->comment);
  gt_str_delete(c->comment_str);
}

static GtStr* comment_node_get_idstr(GtGenomeNode *gn)
{
  GtCommentNode *c;
  gt_assert(gn);
  c = gt_comment_node_cast(gn);
  return c->comment_str;
}

static GtRange comment_node_get_range(GT_UNUSED GtGenomeNode *gn)
{
  GtRange range;
  range.start = 0;
  range.end = 0;
  return range;
}

static int comment_node_accept(GtGenomeNode *gn, GtNodeVisitor *nv,
                               GtError *err)
{
  GtCommentNode *c;
  gt_error_check(err);
  c = gt_comment_node_cast(gn);
  return gt_node_visitor_visit_comment_node(nv, c, err);
}

const GtGenomeNodeClass* gt_comment_node_class()
{
  static const GtGenomeNodeClass *gnc = NULL;
  gt_class_alloc_lock_enter();
  if (!gnc) {
    gnc = gt_genome_node_class_new(sizeof (GtCommentNode),
                                   comment_node_free,
                                   NULL,
                                   comment_node_get_idstr,
                                   comment_node_get_range,
                                   NULL,
                                   NULL,
                                   comment_node_accept);
  }
  gt_class_alloc_lock_leave();
  return gnc;
}

GtGenomeNode* gt_comment_node_new(const char *comment)
{
  GtGenomeNode *gn = gt_genome_node_create(gt_comment_node_class());
  GtCommentNode *c = gt_comment_node_cast(gn);
  gt_assert(comment);
  c->comment = gt_cstr_dup(comment);
  c->comment_str = gt_str_new_cstr("");
  return gn;
}

const char* gt_comment_node_get_comment(const GtCommentNode *c)
{
  gt_assert(c && c->comment);
  return c->comment;
}

GtCommentNode* gt_comment_node_try_cast(GtGenomeNode *gn)
{
  return gt_genome_node_try_cast(gt_comment_node_class(), gn);
}

GtCommentNode* gt_comment_node_cast(GtGenomeNode *gn)
{
  return gt_genome_node_cast(gt_comment_node_class(), gn);
}
