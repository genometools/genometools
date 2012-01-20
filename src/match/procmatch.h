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

#ifndef PROCMATCH_H
#define PROCMATCH_H

#include "core/types_api.h"

typedef struct
{
  bool dbabsolute;
  unsigned long dbseqnum;
  unsigned long dbstartpos,
         dblen;
  const GtUchar *dbsubstring;
  unsigned long querystartpos,
                querylen,
                distance;
  const void *alignment;
} GtIdxMatch;

typedef void (*ProcessIdxMatch)(void *processinfo,
                             const GtIdxMatch *match);

typedef void (*Processresult)(void *,
                              const void *,
                              unsigned long,
                              unsigned long,
                              unsigned long,
                              unsigned long);

#endif
