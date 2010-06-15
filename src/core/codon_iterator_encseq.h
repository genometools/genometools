/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#ifndef CODON_ITERATOR_ENCSEQ_H
#define CODON_ITERATOR_ENCSEQ_H

#include "core/encseq.h"
#include "core/codon_iterator.h"

typedef struct GtCodonIteratorEncseq GtCodonIteratorEncseq;

GtCodonIterator*            gt_codon_iterator_encseq_new(GtEncseq *encseq,
                                                         unsigned long startpos,
                                                         unsigned long length,
                                                         GtError *err);

const GtCodonIteratorClass* gt_codon_iterator_encseq_class(void);
int                         gt_codon_iterator_encseq_unit_test(GtError *err);
#endif
