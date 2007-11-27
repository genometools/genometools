/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef SEQ_H
#define SEQ_H

#include "libgtcore/alpha.h"

typedef struct Seq Seq;

/* stores <seq> pointer */
Seq*          seq_new(const char *seq, unsigned long seqlen, Alpha *seqalpha);
/* takes ownership of <seq> */
Seq*          seq_new_own(char *seq, unsigned long seqlen,
                          Alpha *seqalpha);
/* stores <desc> pointer */
void          seq_set_description(Seq*, const char *desc);
/* takes ownership of <desc> */
void          seq_set_description_own(Seq*, char *desc);
void          seq_set_description(Seq*, const char *desc);
const char*   seq_get_description(Seq*);
const char*   seq_get_orig(const Seq*); /* not '\0' terminated */
const char*   seq_get_encoded(Seq*);
const Alpha*  seq_get_alpha(const Seq*);
unsigned long seq_length(const Seq*);
void          seq_delete(Seq*);

#endif
