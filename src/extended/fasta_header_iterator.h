/*
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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

#ifndef FASTA_HEADER_ITERATOR_H
#define FASTA_HEADER_ITERATOR_H

#include "core/str_array_api.h"
#include "extended/cstr_iterator_rep.h"

/* The <GtFastaHeaderIterator> class, implements <GtCstrIterator> */
typedef struct GtFastaHeaderIterator GtFastaHeaderIterator;

const GtCstrIteratorClass* gt_fasta_header_iterator_class(void);

/* Returns new <GtCstrIterator> which will allow to iterate over the fasta
   headers of file(s) <filenametab>. Returns <NULL> on error. */
GtCstrIterator*            gt_fasta_header_iterator_new(GtStrArray *filenametab,
                                                        GtError *err);

/* Tests if <GSI> is a <GtFastaHeaderIterator> and returns the meta class object
   or <NULL> if the class is not a <GtFastaHeaderIterator>. */
#define                    gt_fasta_header_iterator_cast(GSI) \
  gt_cstr_iterator_cast(gt_fasta_header_iterator_class(), GSI)
#endif
