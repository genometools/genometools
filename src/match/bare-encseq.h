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
#include "core/range_api.h"
#include "core/alphabet_api.h"

typedef struct GtBareSpecialrangeiterator GtBareSpecialrangeiterator;

typedef struct GtBareEncseq GtBareEncseq;

void gt_bare_encseq_delete(GtBareEncseq *bare_encseq);

GtBareEncseq *gt_bare_encseq_new(GtUchar *sequence,GtUword len,
                                 GtUword numofchars);

GtBareEncseq *gt_bare_encseq_parse_new(GtUchar *filecontents,size_t numofbytes,
                                       const GtAlphabet *alphabet,
                                       GtError *err);

const GtUchar *gt_bare_encseq_sequence(const GtBareEncseq *bare_encseq);

GtUword gt_bare_encseq_total_length(const GtBareEncseq *bare_encseq);

GtUword gt_bare_encseq_numofchars(const GtBareEncseq *bare_encseq);

GtUword gt_bare_encseq_charcount(const GtBareEncseq *bare_encseq,
                                 GtUchar idx);

GtUword gt_bare_encseq_specialcharacters(const GtBareEncseq *bare_encseq);

GtUchar gt_bare_encseq_get_encoded_char(const GtBareEncseq *bare_encseq,
                                        GtUword position);

GtBareSpecialrangeiterator* gt_bare_encseq_specialrangeiterator_new(
                       const GtBareEncseq *bare_encseq,
                       bool moveforward);

void gt_bare_encseq_specialrangeiterator_delete(
                       GtBareSpecialrangeiterator *sri);

bool gt_bare_encseq_specialrangeiterator_next(GtBareSpecialrangeiterator *sri,
                                              GtRange *range);

void bare_encseq_convert(GtBareEncseq *bare_encseq,bool forward,bool direct);

#endif
