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

#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "marksubstring.h"

Gtmarksubstring *gt_marksubstring_new(unsigned int numofchars,
                                      unsigned int kmersize,
                                      bool usesuffix,
                                      unsigned int units)
{
  Gtmarksubstring *mark;

  mark = gt_malloc(sizeof (*mark));
  gt_assert(kmersize >= units);
  mark->units = units;
  mark->entries = gt_power_for_small_exponents(numofchars,units);
  if (usesuffix)
  {
    mark->shiftright = 0;
    mark->mask = (GtCodetype) (mark->entries-1);
  } else
  {
    mark->shiftright = GT_MULT2(kmersize - units);
    mark->mask = ~(GtCodetype) 0;
  }
  GT_INITBITTAB(mark->bits,mark->entries);
  mark->size = sizeof (GtBitsequence) * GT_NUMOFINTSFORBITS(mark->entries);
  return mark;
}

void gt_marksubstring_delete(Gtmarksubstring *mark,bool withbits)
{
  if (mark != NULL)
  {
    if (withbits)
    {
      gt_free(mark->bits);
    }
    gt_free(mark);
  }
}

void gt_marksubstring_mark(Gtmarksubstring *mark,GtCodetype code)
{
  code = (code >> (GtCodetype) mark->shiftright) & mark->mask;

  gt_assert(code < mark->entries);
  if (!GT_ISIBITSET(mark->bits,code))
  {
    GT_SETIBIT(mark->bits,code);
  }
}

bool gt_marksubstring_checkmark(const Gtmarksubstring *mark,GtCodetype code)
{
  code = (code >> (GtCodetype) mark->shiftright) & mark->mask;

  return GT_ISIBITSET(mark->bits,code) ? true : false;
}

unsigned int gt_marksubstring_shiftright(const Gtmarksubstring *mark)
{
  return mark->shiftright;
}

unsigned long gt_marksubstring_entries(const Gtmarksubstring *mark)
{
  return mark->entries;
}

size_t gt_marksubstring_size(const Gtmarksubstring *mark)
{
  return mark->size;
}

void gt_marksubstring_bits_null(GT_UNUSED const Gtmarksubstring *mark,bool null)
{
  if (null)
  {
    gt_assert(mark->bits == NULL);
  } else
  {
    gt_assert(mark->bits != NULL);
  }
}

GtBitsequence** gt_marksubstring_bits_address(Gtmarksubstring *mark)
{
  return &mark->bits;
}

void gt_marksubstring_bits_map(Gtmarksubstring *mark,GtBitsequence *bitsmap)
{
  mark->bits = bitsmap;
}
