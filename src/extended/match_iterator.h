/*
  Copyright (c) 2010      Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c)      2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef MATCH_ITERATOR_H
#define MATCH_ITERATOR_H

#include "core/error_api.h"
#include "extended/globalchaining.h"

typedef struct GtMatchIterator GtMatchIterator;
typedef struct GtMatchIteratorClass GtMatchIteratorClass;

typedef enum {
  GT_MATCHER_STATUS_OK,
  GT_MATCHER_STATUS_END,
  GT_MATCHER_STATUS_ERROR
} GtMatchIteratorStatus;

GtMatchIteratorStatus gt_match_iterator_next(GtMatchIterator *mp,
                                             GtFragment *match, GtError *err);
void*                 gt_match_iterator_cast(const GtMatchIteratorClass *mc,
                                             GtMatchIterator *m);
void                   gt_match_iterator_delete(GtMatchIterator*);

#endif
