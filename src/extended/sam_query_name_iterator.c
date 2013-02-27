/*
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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
#include "core/error_api.h"
#include "core/seq_iterator_api.h"
#include "core/unused_api.h"
#include "extended/sam_query_name_iterator.h"
#include "extended/sam_alignment.h"
#include "extended/samfile_iterator.h"
#include "extended/cstr_iterator_rep.h"

struct GtSamQueryNameIterator {
  const GtCstrIterator parent_instance;
  GtSamAlignment *alignment;
  GtSamfileIterator *s_iter;
};

static int gt_sam_query_name_iterator_next(GtCstrIterator *cstr_iterator,
                                      const char **query_name,
                                      GT_UNUSED GtError *err)
{
  int read = 0;
  GtSamQueryNameIterator *sqi = gt_sam_query_name_iterator_cast(cstr_iterator);
  read = gt_samfile_iterator_next(sqi->s_iter, &sqi->alignment);
  while (read > 0 && gt_sam_alignment_is_unmapped(sqi->alignment)) {
    read = gt_samfile_iterator_next(sqi->s_iter, &sqi->alignment);
  }
  if (read > 0) {
    *query_name = gt_sam_alignment_identifier(sqi->alignment);
    return read;
  }
  return 0;
}

static int gt_sam_query_name_iterator_reset(GtCstrIterator *cstr_iterator,
                                       GtError *err)
{
  GtSamQueryNameIterator *sqi = gt_sam_query_name_iterator_cast(cstr_iterator);
  gt_error_check(err);

  return gt_samfile_iterator_reset(sqi->s_iter, err);
}

static void gt_sam_query_name_iterator_delete(GtCstrIterator *cstr_iterator)
{
  GtSamQueryNameIterator *sqi = gt_sam_query_name_iterator_cast(cstr_iterator);
  gt_sam_alignment_delete(sqi->alignment);
}

/* map static local method to interface */
const GtCstrIteratorClass* gt_sam_query_name_iterator_class(void)
{
  static const GtCstrIteratorClass *sic = NULL;
  gt_class_alloc_lock_enter();
  if (!sic) {
    sic = gt_cstr_iterator_class_new(sizeof (GtSamQueryNameIterator),
                                     gt_sam_query_name_iterator_next,
                                     gt_sam_query_name_iterator_reset,
                                     gt_sam_query_name_iterator_delete);
  }
  gt_class_alloc_lock_leave();
  return sic;
}

GtCstrIterator* gt_sam_query_name_iterator_new(GtSamfileIterator *s_iter,
                                        GtError *err)
{
  GtCstrIterator *cstr_iterator =
    gt_cstr_iterator_create(gt_sam_query_name_iterator_class());
  GtSamQueryNameIterator *sqi = gt_sam_query_name_iterator_cast(cstr_iterator);
  sqi->s_iter = s_iter;
  if (gt_sam_query_name_iterator_reset(cstr_iterator, err) != 0)
    return NULL;
  return cstr_iterator;
}
