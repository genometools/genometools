/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef BIOSEQ_ITERATOR_H
#define BIOSEQ_ITERATOR_H

#include "libgtcore/bioseq.h"

typedef struct BioseqIterator BioseqIterator;

/* Create a new BioseqIterator for <seqfile_counter> many <sequence_files>.
   If <seqfile_counter> is 0 use stdin as the only sequence file. */
BioseqIterator* bioseq_iterator_new(int seqfile_counter,
                                    const char **sequence_files);
void            bioseq_iterator_delete(BioseqIterator*);
/* Assign the next <bioseq> if it is available or NULL otherwise.
   Returns -1 in case of failure or 0 otherwise. */
int             bioseq_iterator_next(BioseqIterator*, Bioseq **bioseq, Error*);

#endif
