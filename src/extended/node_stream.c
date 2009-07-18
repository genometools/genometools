/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <stdarg.h>
#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"

struct GtNodeStreamClass {
  size_t size;
  GtNodeStreamFreeFunc free;
  GtNodeStreamNextFunc next;
};

struct GtNodeStreamMembers {
  GtGenomeNode *buffer;
  bool ensure_sorting;
  unsigned int reference_count;
};

const GtNodeStreamClass*
gt_node_stream_class_new(size_t size, GtNodeStreamFreeFunc free,
                         GtNodeStreamNextFunc next)
{
  GtNodeStreamClass *c_class;
  gt_assert(size && next);
  c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free = free;
  c_class->next = next;
  return c_class;
}

GtNodeStream* gt_node_stream_create(const GtNodeStreamClass *nsc,
                                    bool ensure_sorting)
{
  GtNodeStream *ns;
  gt_assert(nsc && nsc->size);
  ns = gt_calloc(1, nsc->size);
  ns->c_class = nsc;
  ns->members = gt_calloc(1, sizeof (GtNodeStreamMembers));
  ns->members->ensure_sorting = ensure_sorting;
  return ns;
}

GtNodeStream* gt_node_stream_ref(GtNodeStream *ns)
{
  gt_assert(ns);
  ns->members->reference_count++;
  return ns;
}

void gt_node_stream_delete(GtNodeStream *ns)
{
  if (!ns) return;
  if (ns->members->reference_count) {
    ns->members->reference_count--;
    return;
  }
  gt_assert(ns->c_class);
  if (ns->c_class->free) ns->c_class->free(ns);
  gt_genome_node_delete(ns->members->buffer);
  gt_free(ns->members);
  gt_free(ns);
}

int gt_node_stream_next(GtNodeStream *ns, GtGenomeNode **gn, GtError *err)
{
  GtGenomeNode *new_node = NULL;
  int had_err = 0;
  gt_assert(ns && ns->c_class && ns->c_class->next);
  gt_error_check(err);
  /* filling */
  if (!ns->members->buffer)
    had_err = ns->c_class->next(ns, &ns->members->buffer, err);
  if (!had_err && ns->members->buffer)
    had_err = ns->c_class->next(ns, &new_node, err);
#ifndef NDEBUG
  /* checking */
  if (!had_err && ns->members->ensure_sorting && ns->members->buffer &&
      new_node) {
    gt_assert(gt_genome_node_compare(&ns->members->buffer, &new_node) <= 0);
  }
#endif
  /* serving */
  if (!had_err) {
    *gn = ns->members->buffer;
    ns->members->buffer = new_node;
  }
  return had_err;
}

int gt_node_stream_pull(GtNodeStream *ns, GtError *err)
{
  GtGenomeNode *gn;
  int had_err;
  gt_assert(ns);
  gt_error_check(err);
  while (!(had_err = gt_node_stream_next(ns, &gn, err)) && gn)
    gt_genome_node_delete(gn);
  return had_err;
}

bool gt_node_stream_is_sorted(GtNodeStream *ns)
{
  gt_assert(ns);
  return ns->members->ensure_sorting;
}

void* gt_node_stream_cast(GT_UNUSED const GtNodeStreamClass *nsc,
                          GtNodeStream *ns)
{
  gt_assert(nsc && ns && ns->c_class == nsc);
  return ns;
}
