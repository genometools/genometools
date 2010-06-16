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

#ifndef CODON_ITERATOR_REP_H
#define CODON_ITERATOR_REP_H

#include <stdio.h>
#include "core/error_api.h"
#include "core/codon_iterator_api.h"

typedef GtCodonIteratorStatus (*CodonIteratorNextFunc)(GtCodonIterator*,
                                                       char*, char*, char*,
                                                       unsigned int*,
                                                       GtError*);
typedef void                  (*CodonIteratorFreeFunc)(GtCodonIterator*);
typedef void                  (*CodonIteratorRewindFunc)(GtCodonIterator*);
typedef unsigned long         (*CodonIteratorCurrentPosFunc)(GtCodonIterator*);
typedef unsigned long         (*CodonIteratorLengthFunc)(GtCodonIterator*);

typedef struct GtCodonIteratorMembers GtCodonIteratorMembers;

struct GtCodonIterator {
  const GtCodonIteratorClass *c_class;
  GtCodonIteratorMembers *pvt;
};

struct GtCodonIteratorMembers {
  unsigned long length,
                curpos,
                startpos;
  unsigned long curframe;
};

const GtCodonIteratorClass* gt_codon_iterator_class_new(size_t size,
                                        CodonIteratorFreeFunc free,
                                        CodonIteratorCurrentPosFunc current_pos,
                                        CodonIteratorLengthFunc length,
                                        CodonIteratorRewindFunc rewind,
                                        CodonIteratorNextFunc next);
GtCodonIterator*            gt_codon_iterator_create(const
                                                         GtCodonIteratorClass*);
void*                       gt_codon_iterator_cast(const GtCodonIteratorClass*,
                                                   GtCodonIterator*);
void*                       gt_codon_iterator_try_cast(const
                                                          GtCodonIteratorClass*,
                                                       GtCodonIterator*);

#endif
