/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef READMODE_DEF_H
#define READMODE_DEF_H
#include "libgtcore/error.h"

typedef enum
{
  Forwardmode = 0,
  Reversemode,
  Complementmode,
  Reversecomplementmode
} Readmode;

#define ISDIRREVERSE(R)    ((R) == Reversemode ||\
                            (R) == Reversecomplementmode)
#define ISDIRCOMPLEMENT(R) ((R) == Complementmode ||\
                            (R) == Reversecomplementmode)

#define COMPLEMENTBASE(B) ((Uchar) 3 - (B))

const char *showreadmode(Readmode readmode);

int parsereadmode(const char *dirargstring,Error *err);

#endif
