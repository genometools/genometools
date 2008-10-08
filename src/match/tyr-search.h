/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef TYR_SEARCH_H
#define TYR_SEARCH_H

#include "core/arraydef.h"
#include "core/str_api.h"
#include "intbits.h"
#include "divmodmul.h"

#define MERBYTES(SL)  (DIV4(SL) + ((MOD4(SL) == 0) ? 0 : 1UL))
#define MERSUFFIX     ".mer"
#define COUNTSSUFFIX  ".mct"
#define BUCKETSUFFIX  ".mbd"
#define EXTRAINTEGERS 2

#define ISBOUNDDEFINED(UDB,IDX)          ISIBITSET(UDB,IDX)
#define SETDEFINEDBOUND(UDB,IDX)         SETIBIT(UDB,IDX)

typedef struct
{
  const void *mappedmbdfileptr;
  const GtStr *indexfilename;
  unsigned int prefixlength;
  unsigned long numofcodes,
                *boundisdefined,
                *bounds;
  uint64_t numberofextractions;
} MBDinfo;

typedef struct
{
  unsigned long idx, value;
} Largecount;

DECLAREARRAYSTRUCT(Largecount);

typedef struct
{
  void *mappedmctfileptr;
  const GtStr *indexfilename;
  Uchar *smallcounts;
  Largecount *largecounts;
  unsigned long numoflargecounts;
} MCTinfo;

typedef struct
{
  void *mappedfileptr;
  const GtStr *indexfilename;
  unsigned int alphasize;
  unsigned long numofmers,
                mersize,
                merbytes;
  Uchar *mertable,
        *lastmer;
} Vmerindex;

typedef struct
{
  Uchar *leftmer,
        *rightmer;
} Merbounds;

#endif
