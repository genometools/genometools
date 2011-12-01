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

#ifndef MARKSUBSTRING_H
#define MARKSUBSTRING_H

#include "core/codetype.h"
#include "core/intbits.h"

typedef struct Gtmarksubstring Gtmarksubstring;

struct Gtmarksubstring
{
  size_t size;
  unsigned int units, shiftright;
  unsigned long entries;
  GtCodetype mask;
  GtBitsequence *bits;
};

#define GT_MARKSUBSTRING_CHECKMARK(MARK,CODE)\
        (marksubstringtmpcode\
           = ((CODE) >> (GtCodetype) (MARK)->shiftright) & (MARK)->mask,\
        GT_ISIBITSET((MARK)->bits,marksubstringtmpcode) ? true : false)

Gtmarksubstring *gt_marksubstring_new(unsigned int numofchars,
                                      unsigned int kmersize, bool usesuffix,
                                      unsigned int units);

void gt_marksubstring_delete(Gtmarksubstring *mark,bool withbits);

void gt_marksubstring_mark(Gtmarksubstring *mark,GtCodetype code);

bool gt_marksubstring_checkmark(const Gtmarksubstring *mark,GtCodetype code);

unsigned int gt_marksubstring_shiftright(const Gtmarksubstring *mark);

unsigned long gt_marksubstring_entries(const Gtmarksubstring *mark);

size_t gt_marksubstring_size(const Gtmarksubstring *mark);

void gt_marksubstring_bits_null(const Gtmarksubstring *mark,bool null);

GtBitsequence** gt_marksubstring_bits_address(Gtmarksubstring *mark);

void gt_marksubstring_bits_map(Gtmarksubstring *mark,
                               GtBitsequence *bitsmap);

#endif
