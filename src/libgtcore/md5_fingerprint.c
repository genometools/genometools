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

#include <ctype.h>
#include "md5.h"
#include "libgtcore/ma.h"
#include "libgtcore/md5_fingerprint.h"
#include "libgtcore/safearith.h"

char *md5_fingerprint(const char *sequence, unsigned long seqlen)
{
  unsigned char output[16];
  char  *upper, *fingerprint;
  unsigned long i;
  /* XXX: this could be done more memory efficient by applying md5 to a reused
     buffer */
  upper = ma_malloc(seqlen * sizeof (char));
  for (i = 0; i < seqlen; i++)
    upper[i] = toupper(sequence[i]);
  md5(upper, safe_cast2long(seqlen), (char*) output);
  ma_free(upper);
  fingerprint = ma_calloc(33, sizeof (char));
  snprintf(fingerprint, 33,
           "%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x",
           output[0], output[1], output[2], output[3], output[4], output[5],
           output[6], output[7], output[8], output[9], output[10], output[11],
           output[12], output[13], output[14], output[15]);
  return fingerprint;
}
