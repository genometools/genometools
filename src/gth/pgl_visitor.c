/*
  Copyright (c) 2008-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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

#include "core/unused_api.h"
#include "gth/pgl_visitor_rep.h"

GthPGLVisitor* gth_pgl_visitor_create(const GthPGLVisitorClass *pglvc)
{
  GthPGLVisitor *pglv;
  gt_assert(pglvc && pglvc->size);
  pglv = gt_calloc(1, pglvc->size);
  pglv->c_class = pglvc;
  return pglv;
}

void* gth_pgl_visitor_cast(GT_UNUSED const GthPGLVisitorClass *pglvc,
                           GthPGLVisitor *pglv)
{
  gt_assert(pglvc && pglv && pglv->c_class == pglvc);
  return pglv;
}

void gth_pgl_visitor_preface(GthPGLVisitor *pglv, unsigned long num_of_pgls)
{
  gt_assert(pglv && pglv->c_class);
  if (pglv->c_class->preface)
    pglv->c_class->preface(pglv, num_of_pgls);
}

void gth_pgl_visitor_set_region_mapping(GthPGLVisitor *pglv,
                                        GtRegionMapping *region_mapping)
{
  gt_assert(pglv && region_mapping && pglv->c_class);
  if (pglv->c_class->set_region_mapping)
    pglv->c_class->set_region_mapping(pglv, region_mapping);
}

void gth_pgl_visitor_visit_pgl(GthPGLVisitor *pglv, GthPGL *pgl,
                               unsigned long pglnum)
{
  gt_assert(pglv && pgl && pglv->c_class && pglv->c_class->visit_pgl);
  pglv->c_class->visit_pgl(pglv, pgl, pglnum);
}

void gth_pgl_visitor_trailer(GthPGLVisitor *pglv)
{
  gt_assert(pglv && pglv->c_class);
  if (pglv->c_class->trailer)
    pglv->c_class->trailer(pglv);
}

void gth_pgl_visitor_delete(GthPGLVisitor *pglv)
{
  if (!pglv) return;
  gt_assert(pglv->c_class);
  if (pglv->c_class->free)
    pglv->c_class->free(pglv);
  gt_free(pglv);
}
