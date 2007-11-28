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

#ifndef FASTABUFFER_H
#define FASTABUFFER_H

#include <stdbool.h>
#include "libgtcore/arraydef.h"
#include "libgtcore/filelengthvalues.h"
#include "libgtcore/queue.h"
#include "libgtcore/strarray.h"
#include "libgtcore/symboldef.h"

typedef struct FastaBuffer FastaBuffer;

FastaBuffer* fastabuffer_new(const StrArray *filenametab,
                             const Uchar *symbolmap, bool plainformat,
                             Filelengthvalues **filelengthtab,
                             Queue *descptr,
                             unsigned long *characterdistribution);
static int   fastabuffer_next(FastaBuffer*, Uchar *val, Error*);
void         fastabuffer_delete(FastaBuffer*);

#include "fastabuffer_imp.h"

#endif
