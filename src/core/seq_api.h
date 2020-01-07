/*
  Copyright (c) 2007-2009 Gordon Gremme <gordon@gremme.org>
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

#ifndef SEQ_API_H
#define SEQ_API_H

#include "core/alphabet_api.h"

/* <GtSeq> is a container for a sequence plus metadata. */
typedef struct GtSeq GtSeq;

/* Create and return a new <GtSeq>, storing the pointer for <seq>
   (of length <seqlen> and with alphabet <seqalpha>). */
GtSeq*            gt_seq_new(const char *seq, GtUword seqlen,
                             GtAlphabet *seqalpha);
/* Like <gt_seq_new()>, but tkes ownership of <seq>. */
GtSeq*            gt_seq_new_own(char *seq, GtUword seqlen,
                                 GtAlphabet *seqalpha);
/* Associates <s> with description <desc>, storing its pointer. */
void              gt_seq_set_description(GtSeq *s, const char *desc);
/* Like <gt_seq_set_description()>, but takes ownership of <desc>. */
void              gt_seq_set_description_own(GtSeq *s, char *desc);
/* Return the description string for <s>. */
const char*       gt_seq_get_description(GtSeq *s);
/* Return the underlying sequence memory for <s>, not guaranteed to be '\0'
   terminated. */
const char*       gt_seq_get_orig(const GtSeq *s);
/* Return the sequence for <s>, encoded using the defined alphabet. */
const GtUchar*    gt_seq_get_encoded(GtSeq *s);
/* Return the alphabet associated with <s>. */
const GtAlphabet* gt_seq_get_alphabet(const GtSeq*);
/* Return the length of the sequence in <s>. */
GtUword           gt_seq_length(const GtSeq *s);
/* Delete <s> and free all associated memory. */
void              gt_seq_delete(GtSeq *s);

#endif
