/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/codon.h"
#include "libgtcore/orf.h"
#include "libgtcore/range.h"
#include "libgtcore/undef.h"

#define START_AMINO 'M'
#define STOP_AMINO  '*'

void determine_ORFs(Array *ranges, unsigned int framenum,
                    const char *frame, unsigned long framelen)
{
  unsigned long i;
  Range orf;
  assert(ranges && framenum < 3 && frame);
  orf.start = UNDEF_ULONG;
  for (i = 0; i < framelen; i++) {
    if (orf.start == UNDEF_ULONG) {
      if (frame[i] == START_AMINO)
        orf.start = i * CODONLENGTH + framenum;
    }
    else {
      if (frame[i] == STOP_AMINO) {
        orf.end = i * CODONLENGTH + framenum + 2;
        array_add(ranges, orf);
        orf.start = UNDEF_ULONG;
      }
    }
  }
}
