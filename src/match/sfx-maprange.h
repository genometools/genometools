/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_MAPRANGE_H
#define SFX_MAPRANGE_H

#include <stdbool.h>
#include <stdlib.h>
#include "core/codetype.h"

#define FROMCODE2SPECIALCODE(CODE,NUMOFCHARS)\
                            (((NUMOFCHARS) == 4U)\
                            ? ((CODE) >> 2)\
                            : (((CODE) - ((NUMOFCHARS)-1)) / (NUMOFCHARS)))

typedef struct
{
  unsigned long mapoffset, mapend;
} GtMappedrange;

void gt_bcktab_mapped_lbrange_get(GtMappedrange *range,
                                  size_t sizeofbasetype,
                                  unsigned long pagesize,
                                  GtCodetype mincode,
                                  GtCodetype maxcode);

unsigned int gt_bcktab_mapped_csrange_get(GtMappedrange *range,
                                          size_t sizeofbasetype,
                                          unsigned long numofallcodes,
                                          unsigned int numofchars,
                                          unsigned long pagesize,
                                          GtCodetype mincode,
                                          GtCodetype maxcode);

#endif
