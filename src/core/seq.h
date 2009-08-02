/*
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/alphabet.h"

typedef struct GtSeq GtSeq;

/* Stores <seq> pointer. */
GtSeq*            gt_seq_new(const char *seq, unsigned long seqlen,
                             GtAlphabet *seqalpha);
/* Takes ownership of <seq>. */
GtSeq*            gt_seq_new_own(char *seq, unsigned long seqlen,
                                 GtAlphabet *seqalpha);
/* Stores <desc> pointer. */
void              gt_seq_set_description(GtSeq*, const char *desc);
/* Takes ownership of <desc>. */
void              gt_seq_set_description_own(GtSeq*, char *desc);
void              gt_seq_set_description(GtSeq*, const char *desc);
const char*       gt_seq_get_description(GtSeq*);
const char*       gt_seq_get_orig(const GtSeq*); /* not '\0' terminated */
const GtUchar*    gt_seq_get_encoded(GtSeq*);
const GtAlphabet* gt_seq_get_alphabet(const GtSeq*);
unsigned long     gt_seq_length(const GtSeq*);
void              gt_seq_delete(GtSeq*);

#endif
