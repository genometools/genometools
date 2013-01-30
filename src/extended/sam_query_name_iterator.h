/*
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg
  Copyright (c) 2010-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef SAM_QUERY_NAME_ITERATOR_H
#define SAM_QUERY_NAME_ITERATOR_H

#include "extended/cstr_iterator_rep.h"
#include "extended/samfile_iterator.h"

/* implements <GtCstrIterator> */
/* iterates over the identifiers of all mapped alignments in a sam/bam file */
typedef struct GtSamQueryNameIterator GtSamQueryNameIterator;

const GtCstrIteratorClass* gt_sam_query_name_iterator_class(void);

/* does not take ownership of <s_iter> */
GtCstrIterator*            gt_sam_query_name_iterator_new(
                                                     GtSamfileIterator *s_iter,
                                                     GtError *err);

#define                    gt_sam_query_name_iterator_cast(GSI) \
  gt_cstr_iterator_cast(gt_sam_query_name_iterator_class(), GSI)
#endif
