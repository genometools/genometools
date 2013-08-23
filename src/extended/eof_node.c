/*
  Copyright (c) 2010 Gordon Gremme <gordon@gremme.org>

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
#include "extended/eof_node_api.h"
#include "extended/genome_node_rep.h"
#include "extended/node_visitor.h"

struct GtEOFNode {
  const GtGenomeNode parent_instance;
};

#define gt_eof_node_cast(eof_node) \
        gt_genome_node_cast(gt_eof_node_class(), eof_node)

static int eof_node_accept(GtGenomeNode *gn, GtNodeVisitor *nv, GtError *err)
{
  GtEOFNode *eof;
  gt_error_check(err);
  eof = gt_eof_node_cast(gn);
  return gt_node_visitor_visit_eof_node(nv, eof, err);
}

const GtGenomeNodeClass* gt_eof_node_class()
{
  static const GtGenomeNodeClass *gnc = NULL;
  gt_class_alloc_lock_enter();
  if (!gnc) {
    gnc = gt_genome_node_class_new(sizeof (GtEOFNode),
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL,
                                   eof_node_accept);
  }
  gt_class_alloc_lock_leave();
  return gnc;
}

GtGenomeNode* gt_eof_node_new(void)
{
  return gt_genome_node_create(gt_eof_node_class());
}

GtEOFNode* gt_eof_node_try_cast(GtGenomeNode *gn)
{
  GtEOFNode *en = gt_genome_node_try_cast(gt_eof_node_class(), gn);
  return en;
}
