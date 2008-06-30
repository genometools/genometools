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

#include "libgtcore/ma.h"
#include "libgtcore/mathsupport.h"
#include "libgtext/shredder.h"

struct Shredder {
  Bioseq *bioseq;
  unsigned long minlength,
                maxlength,
                overlap,
                seqnum,
                pos;
};

Shredder* shredder_new(Bioseq *bioseq, unsigned long minlength,
                                       unsigned long maxlength)
{
  Shredder *shredder = ma_calloc(1, sizeof *shredder);
  assert(bioseq && minlength && minlength <= maxlength);
  shredder->bioseq = bioseq;
  shredder->minlength = minlength;
  shredder->maxlength = maxlength;
  return shredder;
}

void shredder_delete(Shredder *shredder)
{
  if (!shredder) return;
  ma_free(shredder);
}

void shredder_set_overlap(Shredder *shredder, unsigned long overlap)
{
  assert(shredder);
  shredder->overlap = overlap;
}

const char* shredder_shred(Shredder *shredder, unsigned long *fragment_length,
                           Str *desc)
{
  assert(shredder && fragment_length);
  if (shredder->seqnum < bioseq_number_of_sequences(shredder->bioseq)) {
    unsigned long seqlen, fraglen;
    const char *frag;
    seqlen = bioseq_get_sequence_length(shredder->bioseq, shredder->seqnum);
    fraglen = (shredder->maxlength == shredder->minlength
               ? 0 : rand_max(shredder->maxlength - shredder->minlength))
              + shredder->minlength;
    assert(fraglen >= shredder->minlength);
    frag = bioseq_get_sequence(shredder->bioseq, shredder->seqnum)
           + shredder->pos;
    if (shredder->pos + fraglen > seqlen)
      fraglen = seqlen - shredder->pos;
    *fragment_length = fraglen;
    str_append_cstr(desc, bioseq_get_description(shredder->bioseq,
                                                 shredder->seqnum));
    assert(shredder->pos + fraglen <= seqlen);
    if (shredder->pos + fraglen == seqlen) { /* last fragment */
      shredder->seqnum++;
      shredder->pos = 0;
    }
    else {
      if (fraglen > shredder->overlap)
        shredder->pos += fraglen - shredder->overlap;
      else
        shredder->pos++; /* go at least one base further each step */
    }
    return frag;
  }
  return NULL;
}
