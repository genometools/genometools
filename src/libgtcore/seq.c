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

#include <assert.h>
#include "libgtcore/ma.h"
#include "libgtcore/seq.h"
#include "libgtcore/xansi.h"

struct Seq {
  char *seq, *description;
  char *encoded_seq;
  unsigned long seqlen;
  bool own_seq,
       own_description;
  Alpha *seqalpha;
};

Seq* seq_new(const char *seq, unsigned long seqlen, Alpha *seqalpha)
{
  Seq *s;
  assert(seq && seqalpha);
  s = ma_calloc(1, sizeof (Seq));
  s->seq = (char*) seq;
  s->seqlen = seqlen;
  s->seqalpha = alpha_ref(seqalpha);
  return s;
}

Seq* seq_new_own(char* seq, unsigned long seqlen, Alpha *seqalpha)
{
  Seq *s = seq_new(seq, seqlen, seqalpha);
  s->own_seq = true;
  return s;
}

void seq_set_description(Seq *s, const char *desc)
{
  assert(s);
  s->description = (char*) desc;
}

void seq_set_description_own(Seq *s, char *desc)
{
  assert(s);
  if (s->description && s->own_description)
    ma_free(s->description);
  s->description = desc;
  s->own_description = true;
}

const char* seq_get_description(Seq *s)
{
  assert(s);
  return s->description;
}

const char* seq_get_orig(const Seq *s)
{
  assert(s);
  return s->seq;
}

const char* seq_get_encoded(Seq *s)
{
  assert(s);
  if (!s->encoded_seq) {
    s->encoded_seq = ma_malloc(sizeof (char) * (s->seqlen+1));
    alpha_encode_seq(s->seqalpha, s->encoded_seq, (char*) s->seq, s->seqlen);
    s->encoded_seq[s->seqlen] = '\0';
  }
  return s->encoded_seq;
}

const Alpha* seq_get_alpha(const Seq *s)
{
  assert(s);
  return s->seqalpha;
}

unsigned long seq_length(const Seq *s)
{
  assert(s);
  return s->seqlen;
}

void seq_delete(Seq *s)
{
  if (!s) return;
  if (s->own_seq)
    ma_free(s->seq);
  if (s->own_description)
    ma_free(s->description);
  ma_free(s->encoded_seq);
  alpha_delete(s->seqalpha);
  ma_free(s);
}
