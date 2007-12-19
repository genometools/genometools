/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#ifndef DATAALIGN_H
#define DATAALIGN_H

#include <stdlib.h>
#include "minmax.h"

/*
 * functionality to layout data at correct alignment (=> fewer
 * individual mallocs)
 */
enum {
  MAX_ALIGN_REQUIREMENT = 8,
  MIN_ALIGN_REQUREMENT = 4,
};

static inline unsigned long long
roundUp(unsigned long long v, unsigned long multipleOf)
{
  return v - v%multipleOf + multipleOf * (v%multipleOf?1:0);
}

static inline size_t
offsetAlign(size_t offset, size_t sizeOfVal2Align)
{
  size_t alignBase = MAX(MIN_ALIGN_REQUREMENT,
                         MIN(sizeOfVal2Align, MAX_ALIGN_REQUIREMENT));
  return roundUp(offset, alignBase);
}

#endif
