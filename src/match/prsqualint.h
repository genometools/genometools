/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef PRSQUALINT_H
#define PRSQUALINT_H

#include "core/error_api.h"

#define BESTCHARACTER        'b'
#define PERCENTAWAYCHARACTER 'p'

typedef enum
{
  Qualabsolute,
  Qualpercentaway,
  Qualbestof
} Qualificationtag;

typedef struct
{
  Qualificationtag qualtag;
  unsigned long integervalue;
} Qualifiedinteger;

Qualifiedinteger *gt_parsequalifiedinteger(const char *option,
                                           const char *lparam,
                                           GtError *err);

#endif
