/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#ifndef EXTEND_OFFSET_H
#define EXTEND_OFFSET_H

#include <stdbool.h>
#include "core/readmode_api.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/assert_api.h"
#include "core/readmode.h"

inline bool gt_extend_read_seq_left2right(bool rightextension,
                                          GtReadmode readmode)
{
  return (rightextension && !GT_ISDIRREVERSE(readmode)) ||
         (!rightextension && GT_ISDIRREVERSE(readmode)) ? true : false;
}

inline GtUword gt_extend_offset(bool rightextension,
                                GtReadmode readmode,
                                GtUword totallength,
                                GtUword startpos,
                                GtUword len,
                                GT_UNUSED GtUword totallength_undef)
{
  GtUword offset;

  if (rightextension)
  {
    if (GT_ISDIRREVERSE(readmode))
    {
      gt_assert(totallength != totallength_undef && startpos < totallength);
      offset = totallength - 1 - startpos;
    } else
    {
      offset = startpos;
    }
  } else
  {
    if (GT_ISDIRREVERSE(readmode))
    {
      gt_assert(totallength != totallength_undef &&
                startpos + totallength >= len);
      offset = startpos + totallength - len;
    } else
    {
      gt_assert(startpos + len > 0);
      offset = startpos + len - 1;
    }
  }
  return offset;
}
#endif
