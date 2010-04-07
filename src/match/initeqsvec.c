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

#include <limits.h>
#include "core/assert_api.h"
#include "core/chardef.h"
#include "initeqsvec.h"

void gt_initeqsvector(unsigned long *eqsvector,
                   unsigned long eqslen,
                   const GtUchar *pattern,
                   unsigned long patternlength)
{
  unsigned long *eptr, shiftmask;
  const GtUchar *pptr;

  gt_assert(patternlength <=
            (unsigned long) (CHAR_BIT * sizeof (unsigned long)));
  for (eptr = eqsvector; eptr < eqsvector + eqslen; eptr++)
  {
    *eptr = 0;
  }
  for (pptr = pattern, shiftmask = 1UL;
       pptr < pattern + patternlength && shiftmask != 0;
       pptr++, shiftmask <<= 1)
  {
    gt_assert(*pptr != (GtUchar) SEPARATOR);
    if (*pptr != (GtUchar) WILDCARD)
    {
      eqsvector[(unsigned long) *pptr] |= shiftmask;
    }
  }
}

void gt_initeqsvectorrev(unsigned long *eqsvectorrev,
                      unsigned long eqslen,
                      const GtUchar *pattern,
                      unsigned long patternlength)
{
  unsigned long *eptr, shiftmask;
  const GtUchar *pptr;

  gt_assert(patternlength <=
            (unsigned long) (CHAR_BIT * sizeof (unsigned long)));
  for (eptr = eqsvectorrev; eptr < eqsvectorrev + eqslen; eptr++)
  {
    *eptr = 0;
  }
  for (pptr = pattern+patternlength-1, shiftmask = 1UL;
       pptr >= pattern && shiftmask != 0;
       pptr--, shiftmask <<= 1)
  {
    gt_assert(*pptr != (GtUchar) SEPARATOR);
    if (*pptr != (GtUchar) WILDCARD)
    {
      eqsvectorrev[(unsigned long) *pptr] |= shiftmask;
    }
  }
}
