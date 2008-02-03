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

#ifndef SUBSTRITER_H
#define SUBSTRITER_H

#include <inttypes.h>
#include <stdbool.h>
#include "libgtcore/error.h"
#include "libgtcore/strarray.h"
#include "libgtcore/symboldef.h"
#include "alphadef.h"

typedef struct Substriter Substriter;

typedef struct
{
  const Uchar *querystart, *queryptr;
  unsigned int currentcode;
  bool definedcurrentcode;
  unsigned long remaining;
  char *desc;
} Substring;

Substriter *substriter_new(const StrArray *queryfilenames,
                           const Alphabet *alphabet,
                           unsigned int qvalue);

int substriter_next(Substring *substring,Substriter *substriter,Error *err);

uint64_t substriter_unitnum(const Substriter *substriter);

void substriter_delete(Substriter **substriter);

#endif
