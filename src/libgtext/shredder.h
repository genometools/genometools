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

#ifndef SHREDDER_H
#define SHREDDER_H

#include "libgtcore/bioseq.h"

typedef struct Shredder Shredder;

/* Create new Shredder for sequences in <bioseq>. The produced fragments will
   have at least length <minlength> and at most length <maxlength>. */
Shredder*   shredder_new(Bioseq *bioseq, unsigned long minlength,
                                         unsigned long maxlength);
void        shredder_delete(Shredder*);
/* Set the <overlap> between shredded fragments, the default is 0. */
void        shredder_set_overlap(Shredder*, unsigned long overlap);
/* Return the next shredded fragment or NULL if no additional fragment is
   available. The length of the fragment is stored in <fragment_length> and the
   description of the corresponding sequence is appended to <desc>. */
const char* shredder_shred(Shredder*, unsigned long *fragment_length,
                           Str *desc);

#endif
