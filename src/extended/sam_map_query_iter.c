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
#include "core/error_api.h"
#include "core/seqiterator.h"
#include "core/unused_api.h"
#include "extended/sam_map_query_iter.h"
#include "extended/sam_alignment.h"
#include "extended/samfile_iterator.h"
#include "extended/string_iter.h"
#include "extended/string_iter_rep.h"

struct GtSamMapQueryIter {
  const GtStringIter parent_instance;
  GtSamAlignment *alignment;
  GtSamfileIterator *s_iter;
};

static int gt_sam_map_query_iter_next(GtStringIter *str_iter,
                                      const char **query_name,
                                      GT_UNUSED GtError *err)
{
  int read = 0;
  GtSamMapQueryIter *sqi = gt_sam_map_query_iter_cast(str_iter);
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

static int gt_sam_map_query_iter_reset(GtStringIter *str_iter,
                                       GtError *err)
{
  GtSamMapQueryIter *sqi = gt_sam_map_query_iter_cast(str_iter);
  gt_error_check(err);

  return gt_samfile_iterator_reset(sqi->s_iter, err);
}

static void gt_sam_map_query_iter_delete(GtStringIter *str_iter)
{
  GtSamMapQueryIter *sqi = gt_sam_map_query_iter_cast(str_iter);
  gt_sam_alignment_delete(sqi->alignment);
}

/* map static local method to interface */
const GtStringIterClass* gt_sam_map_query_iter_class(void)
{
  static const GtStringIterClass *sic = NULL;
  if (sic == NULL) {
    sic = gt_string_iter_class_new(sizeof (GtSamMapQueryIter),
                                   gt_sam_map_query_iter_next,
                                   gt_sam_map_query_iter_reset,
                                   gt_sam_map_query_iter_delete);
  }
  return sic;
}

GtStringIter* gt_sam_map_query_iter_new(GtSamfileIterator *s_iter,
                                        GtError *err)
{
  GtStringIter *str_iter = gt_string_iter_create(gt_sam_map_query_iter_class());
  GtSamMapQueryIter *sqi = gt_sam_map_query_iter_cast(str_iter);
  sqi->s_iter = s_iter;
  if (gt_sam_map_query_iter_reset(str_iter, err) != 0)
    return NULL;
  return str_iter;
}
