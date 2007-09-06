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

#include <stdio.h>
#include <stdlib.h>
#include <libgtcore/safeop.h>

long safe_cast_to_long_type(unsigned long ulong, const char *filename, int line)
{
  if (ulong > ~(1UL << ((8UL * sizeof (unsigned long)) - 1))) {
    fprintf(stderr, "%s, l.%d: illegal cast to long (ulong=%lu)\n",
            filename, line, ulong);
    exit(EXIT_FAILURE);
  }
  return (long) ulong;
}

unsigned long safe_cast_to_ulong_type(long slong, const char *filename,
                                      int line)
{
  if (slong < 0) {
    fprintf(stderr, "%s, l.%d: illegal cast to unsigned long (slong=%ld)\n",
            filename, line, slong);
    exit(EXIT_FAILURE);
  }
  return (unsigned long) slong;
}

