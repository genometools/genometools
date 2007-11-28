/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQITERATOR_H
#define SEQITERATOR_H

#include <inttypes.h>
#include "libgtcore/queue.h"
#include "libgtcore/strarray.h"
#include "libgtcore/symboldef.h"

typedef struct SeqIterator SeqIterator;

SeqIterator* seqiterator_new(const StrArray *filenametab,
                             const Uchar *symbolmap, bool withsequence);
int          seqiterator_next(SeqIterator*, const Uchar **sequence,
                              unsigned long *len, char **desc, Error*);
const unsigned long long*
             seqiterator_getcurrentcounter(SeqIterator*, unsigned long long);
void         seqiterator_delete(SeqIterator*);

#endif
