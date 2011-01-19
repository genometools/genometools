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

#include "core/ma.h"
#include "core/seq.h"

struct GtSeq {
  char *seq, *description;
  GtUchar *encoded_seq;
  unsigned long seqlen;
  bool own_seq,
       own_description;
  GtAlphabet *seqalpha;
};

GtSeq* gt_seq_new(const char *seq, unsigned long seqlen, GtAlphabet *seqalpha)
{
  GtSeq *s;
  gt_assert(seq && seqalpha);
  s = gt_calloc(1, sizeof (GtSeq));
  s->seq = (char*) seq;
  s->seqlen = seqlen;
  s->seqalpha = gt_alphabet_ref(seqalpha);
  return s;
}

GtSeq* gt_seq_new_own(char* seq, unsigned long seqlen, GtAlphabet *seqalpha)
{
  GtSeq *s = gt_seq_new(seq, seqlen, seqalpha);
  s->own_seq = true;
  return s;
}

void gt_seq_set_description(GtSeq *s, const char *desc)
{
  gt_assert(s);
  s->description = (char*) desc;
}

void gt_seq_set_description_own(GtSeq *s, char *desc)
{
  gt_assert(s);
  if (s->description && s->own_description)
    gt_free(s->description);
  s->description = desc;
  s->own_description = true;
}

const char* gt_seq_get_description(GtSeq *s)
{
  gt_assert(s);
  return s->description;
}

const char* gt_seq_get_orig(const GtSeq *s)
{
  gt_assert(s);
  return s->seq;
}

const GtUchar* gt_seq_get_encoded(GtSeq *s)
{
  gt_assert(s);
  if (!s->encoded_seq) {
    s->encoded_seq = gt_malloc(sizeof (char) * (s->seqlen+1));
    gt_alphabet_encode_seq(s->seqalpha, s->encoded_seq, (char*) s->seq,
                           s->seqlen);
    s->encoded_seq[s->seqlen] = '\0';
  }
  return s->encoded_seq;
}

const GtAlphabet* gt_seq_get_alphabet(const GtSeq *s)
{
  gt_assert(s);
  return s->seqalpha;
}

unsigned long gt_seq_length(const GtSeq *s)
{
  gt_assert(s);
  return s->seqlen;
}

void gt_seq_delete(GtSeq *s)
{
  if (!s) return;
  if (s->own_seq)
    gt_free(s->seq);
  if (s->own_description)
    gt_free(s->description);
  gt_free(s->encoded_seq);
  gt_alphabet_delete(s->seqalpha);
  gt_free(s);
}
