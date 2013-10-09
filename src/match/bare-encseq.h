/*
  Copyright (c) 2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#ifndef BARE_ENCSEQ_H
#define BARE_ENCSEQ_H

#include "core/error_api.h"
#include "core/types_api.h"

typedef struct
{
  GtUword start,
          length;
} GtSainSpecialrange;

GT_DECLAREARRAYSTRUCT(GtSainSpecialrange);

typedef struct GtBareEncseq GtBareEncseq;

void gt_bare_encseq_delete(GtBareEncseq *bare_encseq);

GtBareEncseq *gt_bare_encseq_new(GtUchar *filecontents,size_t numofbytes,
                                 GtError *err);

const GtArrayGtSainSpecialrange *gt_bare_encseq_specialranges(
                                           const GtBareEncseq *bare_encseq);

const GtUchar *gt_bare_encseq_sequence(const GtBareEncseq *bare_encseq);

GtUword gt_bare_encseq_total_length(const GtBareEncseq *bare_encseq);

GtUword gt_bare_encseq_numofchars(const GtBareEncseq *bare_encseq);

GtUword gt_bare_encseq_charcount(const GtBareEncseq *bare_encseq,
                                 GtUchar idx);

GtUword gt_bare_encseq_specialcharacters(const GtBareEncseq *bare_encseq);

#endif
