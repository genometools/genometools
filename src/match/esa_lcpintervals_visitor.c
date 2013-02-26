/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include "core/class_alloc_lock.h"
#include "core/unused_api.h"
#include "esa_visitor_rep.h"
#include "esa_lcpintervals_visitor.h"

struct GtESALcpintervalsVisitor {
  const GtESAVisitor parent_instance;
};

#define gt_esa_lcpitvs_visitor_cast(GV)\
        gt_esa_visitor_cast(gt_esa_lcpitvs_visitor_class(), GV)

static int gt_esa_lcpitvs_visitor_processleafedge(GT_UNUSED GtESAVisitor *ev,
                                                  bool firstsucc,
                                                  unsigned long fd,
                                                  GT_UNUSED unsigned long flb,
                                                  GT_UNUSED
                                                         GtESAVisitorInfo *info,
                                                  unsigned long leafnumber,
                                                  GT_UNUSED GtError *err)

{
  printf("L %c %lu %lu %lu\n", firstsucc ? '1' : '0', fd, flb, leafnumber);
  return 0;
}

static int gt_esa_lcpitvs_visitor_processbranchingedge(
                                                    GT_UNUSED GtESAVisitor *ev,
                                                    bool firstsucc,
                                                    unsigned long fd,
                                                    unsigned long flb,
                                                    GT_UNUSED
                                                        GtESAVisitorInfo *finfo,
                                                    unsigned long sd,
                                                    unsigned long slb,
                                                    GT_UNUSED unsigned long srb,
                                                    GT_UNUSED
                                                        GtESAVisitorInfo *sinfo,
                                                    GT_UNUSED GtError *err)
{
  printf("B %c %lu %lu %lu %lu\n", firstsucc ? '1' : '0', fd, flb, sd, slb);
  return 0;
}

static const GtESAVisitorClass* gt_esa_lcpitvs_visitor_class(void)
{
  static const GtESAVisitorClass *esc = NULL;
  gt_class_alloc_lock_enter();
  if (!esc) {
    esc = gt_esa_visitor_class_new(sizeof (GtESALcpintervalsVisitor),
                                   NULL,
                                   gt_esa_lcpitvs_visitor_processleafedge,
                                   gt_esa_lcpitvs_visitor_processbranchingedge,
                                   NULL,
                                   NULL,
                                   NULL);
  }
  gt_class_alloc_lock_leave();
  return esc;
}

GtESAVisitor* gt_esa_lcpitvs_visitor_new(void)
{
  GtESAVisitor *ev = gt_esa_visitor_create(gt_esa_lcpitvs_visitor_class());
  return ev;
}
