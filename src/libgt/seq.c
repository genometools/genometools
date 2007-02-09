/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "seq.h"
#include "xansi.h"

struct Seq {
  const char *seq, *description;
  char *encoded_seq;
  unsigned long seqlen;
  Alpha *seqalpha;
};

Seq* seq_new(const char *seq, unsigned long seqlen, Alpha *seqalpha)
{
  Seq *s;
  assert(seq && seqalpha);
  s = xcalloc(1, sizeof(Seq));
  s->seq = seq;
  s->seqlen = seqlen;
  s->seqalpha = alpha_ref(seqalpha);
  return s;
}

void seq_set_description(Seq *s, const char *desc)
{
  assert(s);
  s->description = desc;
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
    s->encoded_seq = xmalloc(sizeof(char) * (s->seqlen+1));
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

void seq_free(Seq *s)
{
  if (!s) return;
  free(s->encoded_seq);
  alpha_free(s->seqalpha);
  free(s);
}
