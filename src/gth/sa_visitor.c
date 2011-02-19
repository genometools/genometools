/*
  Copyright (c) 2008-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "gth/sa_visitor_rep.h"

GthSAVisitor* gth_sa_visitor_create(const GthSAVisitorClass *savc)
{
  GthSAVisitor *sav;
  gt_assert(savc && savc->size);
  sav = gt_calloc(1, savc->size);
  sav->c_class = savc;
  return sav;
}

void* gth_sa_visitor_cast(GT_UNUSED const GthSAVisitorClass *savc,
                          GthSAVisitor *sav)
{
  gt_assert(savc && sav && sav->c_class == savc);
  return sav;
}

void gth_sa_visitor_preface(GthSAVisitor *sav)
{
  gt_assert(sav && sav->c_class);
  if (sav->c_class->preface)
    sav->c_class->preface(sav);
}

void gth_sa_visitor_visit_sa(GthSAVisitor *sav, GthSA *sa)
{
  gt_assert(sav && sa && sav->c_class && sav->c_class->visit_sa);
  sav->c_class->visit_sa(sav, sa);
}

void gth_sa_visitor_trailer(GthSAVisitor *sav, unsigned long num_of_sas)
{
  gt_assert(sav && sav->c_class);
  if (sav->c_class->trailer)
    sav->c_class->trailer(sav, num_of_sas);
}

void gth_sa_visitor_delete(GthSAVisitor *sav)
{
  if (!sav) return;
  gt_assert(sav->c_class);
  if (sav->c_class->free)
    sav->c_class->free(sav);
  gt_free(sav);
}
