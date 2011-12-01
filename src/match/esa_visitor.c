/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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
#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "match/esa_visitor_rep.h"

/* the <GtESAVisitor> interface */
struct GtESAVisitorClass {
  size_t size;
  GtESAVisitorFreeFunc free;
  GtESAVisitorLeafEdgeFunc leafedge_func;
  GtESAVisitorBranchingEdgeFunc branchingedge_func;
  GtESAVisitorLCPIntervalFunc lcpinterval_func;
  GtESAVisitorInfoCreatorFunc create_info_func;
  GtESAVisitorInfoDestructorFunc delete_info_func;
};

const GtESAVisitorClass*
gt_esa_visitor_class_new(size_t size,
                         GtESAVisitorFreeFunc free,
                         GtESAVisitorLeafEdgeFunc leafedge_func,
                         GtESAVisitorBranchingEdgeFunc branchingedge_func,
                         GtESAVisitorLCPIntervalFunc lcpinterval_func,
                         GtESAVisitorInfoCreatorFunc create_info_func,
                         GtESAVisitorInfoDestructorFunc delete_info_func)
{
  GtESAVisitorClass *c_class;
  gt_assert(size);
  c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free = free;
  c_class->leafedge_func = leafedge_func;
  c_class->branchingedge_func = branchingedge_func;
  c_class->lcpinterval_func = lcpinterval_func;
  c_class->create_info_func = create_info_func;
  c_class->delete_info_func = delete_info_func;
  return c_class;
}

GtESAVisitor* gt_esa_visitor_create(const GtESAVisitorClass *evc)
{
  GtESAVisitor *ev;
  gt_assert(evc && evc->size);
  ev = gt_calloc((size_t) 1, evc->size);
  ev->c_class = evc;
  return ev;
}

void* gt_esa_visitor_cast(GT_UNUSED const GtESAVisitorClass *evc,
                           GtESAVisitor *ev)
{
  gt_assert(evc && ev && ev->c_class == evc);
  return ev;
}

GtESAVisitorInfo* gt_esa_visitor_info_new(GtESAVisitor *ev)
{
  gt_assert(ev && ev->c_class);
  if (ev->c_class->create_info_func != NULL)
    return ev->c_class->create_info_func(ev);
  return NULL;
}

void gt_esa_visitor_info_delete(GtESAVisitorInfo *info, GtESAVisitor *ev)
{
  gt_assert(ev && ev->c_class);
  if (ev->c_class->delete_info_func != NULL)
    ev->c_class->delete_info_func(info, ev);
}

int gt_esa_visitor_visit_leaf_edge(GtESAVisitor *ev,
                                   bool firstsucc,
                                   unsigned long fd,
                                   unsigned long flb,
                                   GtESAVisitorInfo *info,
                                   unsigned long leafnumber,
                                   GtError *err)
{
  gt_error_check(err);
  gt_assert(ev && ev->c_class);
  if (ev->c_class->leafedge_func != NULL)
    return ev->c_class->leafedge_func(ev, firstsucc, fd, flb, info, leafnumber,
                                      err);
  return 0;
}

int gt_esa_visitor_visit_branching_edge(GtESAVisitor *ev,
                                        bool firstsucc,
                                        unsigned long fd,
                                        unsigned long flb,
                                        GtESAVisitorInfo *finfo,
                                        unsigned long sd,
                                        unsigned long slb,
                                        unsigned long srb,
                                        GtESAVisitorInfo *sinfo,
                                        GtError *err)
{
  gt_error_check(err);
  gt_assert(ev && ev->c_class);
  if (ev->c_class->branchingedge_func != NULL)
    return ev->c_class->branchingedge_func(ev, firstsucc, fd, flb, finfo, sd,
                                           slb, srb, sinfo, err);
  return 0;
}

int gt_esa_visitor_visit_lcp_interval(GtESAVisitor *ev,
                                      unsigned long lcp,
                                      unsigned long lb,
                                      unsigned long rb,
                                      GtESAVisitorInfo *info,
                                      GtError *err)
{
  gt_error_check(err);
  gt_assert(ev && ev->c_class);
  if (ev->c_class->lcpinterval_func != NULL)
    return ev->c_class->lcpinterval_func(ev, lcp, lb, rb, info, err);
  return 0;
}

void gt_esa_visitor_delete(GtESAVisitor *ev)
{
  if (!ev) return;
  gt_assert(ev->c_class);
  if (ev->c_class->free != NULL)
    ev->c_class->free(ev);
  gt_free(ev);
}
