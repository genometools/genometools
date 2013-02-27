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
#include <string.h>
#include "md5.h"
#include "core/ma.h"
#include "core/md5_encoder_api.h"
#include "core/md5_fingerprint_api.h"
#include "core/safearith.h"

char* gt_md5_fingerprint(const char *sequence, unsigned long seqlen)
{
  unsigned char output[16];
  char buf[64];
  char *fingerprint;
  GtMD5Encoder *enc;
  unsigned long i, pos = 0;

  enc = gt_md5_encoder_new();
  for (i = 0; i < seqlen; i++) {
    if (pos == 64) {
      gt_md5_encoder_add_block(enc, buf, 64);
      pos = 0;
    }
    buf[pos++] = toupper(sequence[i]);
  }
  gt_md5_encoder_add_block(enc, buf, pos);
  fingerprint = gt_calloc(33, sizeof (char));
  gt_md5_encoder_finish(enc, output, fingerprint);

  gt_md5_encoder_delete(enc);
  return fingerprint;
}
