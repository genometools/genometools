/*
  Copyright (c) 2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "extended/meta_node_api.h"
#include "extended/genome_node_rep.h"

struct GtMetaNode
{
  const GtGenomeNode parent_instance;
  char *meta_directive,
       *meta_data;
  GtStr *meta_str; /* used in gt_meta_node_get_idstr() */
};

static void meta_node_free(GtGenomeNode *gn)
{
  GtMetaNode *m = gt_meta_node_cast(gn);
  gt_assert(m && m->meta_directive && m->meta_data);
  gt_free(m->meta_directive);
  gt_free(m->meta_data);
  gt_str_delete(m->meta_str);
}

static GtStr* meta_node_get_idstr(GtGenomeNode *gn)
{
  GtMetaNode *m;
  gt_assert(gn);
  m = gt_meta_node_cast(gn);
  return m->meta_str;
}

static GtRange meta_node_get_range(GT_UNUSED GtGenomeNode *gn)
{
  GtRange range;
  range.start = 0;
  range.end = 0;
  return range;
}

static int meta_node_accept(GtGenomeNode *gn, GtNodeVisitor *nv, GtError *err)
{
  GtMetaNode *m;
  gt_error_check(err);
  m = gt_meta_node_cast(gn);
  return gt_node_visitor_visit_meta_node(nv, m, err);
}

const GtGenomeNodeClass* gt_meta_node_class()
{
  static const GtGenomeNodeClass *gnc = NULL;
  gt_class_alloc_lock_enter();
  if (!gnc) {
    gnc = gt_genome_node_class_new(sizeof (GtMetaNode),
                                   meta_node_free,
                                   NULL,
                                   meta_node_get_idstr,
                                   meta_node_get_range,
                                   NULL,
                                   NULL,
                                   meta_node_accept);
  }
  gt_class_alloc_lock_leave();
  return gnc;
}

GtGenomeNode* gt_meta_node_new(const char *meta_directive,
                               const char *meta_data)
{
  GtGenomeNode *gn = gt_genome_node_create(gt_meta_node_class());
  GtMetaNode *m = gt_meta_node_cast(gn);
  gt_assert(meta_directive && meta_data);
  m->meta_directive = gt_cstr_dup(meta_directive);
  m->meta_data = gt_cstr_dup(meta_data);
  m->meta_str = gt_str_new_cstr("");
  return gn;
}

const char* gt_meta_node_get_directive(const GtMetaNode *m)
{
  gt_assert(m && m->meta_directive);
  return m->meta_directive;
}

const char* gt_meta_node_get_data(const GtMetaNode *m)
{
  gt_assert(m && m->meta_data);
  return m->meta_data;
}
