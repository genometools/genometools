/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/ma_api.h"
#include "core/intbits.h"
#include "sfx-suffixgetset.h"
#include "radixsort_str.h"
#include "sfx-radixsort.h"
#include "sfx-lcpvalues.h"

void gt_sfx_radixsort_str(GtRadixsortstringinfo *rsi,
                          unsigned long depth,
                          unsigned int sortmaxdepth,
                          unsigned long subbucketleft,
                          unsigned long width,
                          GtSuffixsortspace *sssp,
                          GtLcpvalues *lcpvalues)
{
  unsigned long idx, *suffixes;
  GtSuffixsortspace_exportptr *exportptr
    = gt_suffixsortspace_exportptr(subbucketleft, sssp);
  bool allocated = false;

  if (exportptr->ulongtabsectionptr != NULL)
  {
    suffixes = exportptr->ulongtabsectionptr;
  } else
  {
    suffixes = gt_malloc(sizeof (*suffixes) * width);
    allocated = true;
    for (idx = 0; idx < width; idx++)
    {
      suffixes[idx] = (unsigned long) exportptr->uinttabsectionptr[idx];
    }
  }
  gt_radixsort_str_eqlen(rsi,
                         suffixes,
                         lcpvalues,
                         subbucketleft,
                         depth,
                         (unsigned long) sortmaxdepth,
                         width);
  if (allocated)
  {
    gt_assert(exportptr->uinttabsectionptr != NULL);
    for (idx = 0; idx < width; idx++)
    {
      exportptr->uinttabsectionptr[idx] = (uint32_t) suffixes[idx];
      if (exportptr->uinttabsectionptr[idx] == 0)
      {
        gt_suffixsortspace_updatelongest(sssp,subbucketleft + idx);
      }
    }
    gt_free(suffixes);
  } else
  {
    gt_assert(exportptr->ulongtabsectionptr != NULL);
    for (idx = 0; idx < width; idx++)
    {
      if (exportptr->ulongtabsectionptr[idx] == 0)
      {
        gt_suffixsortspace_updatelongest(sssp,subbucketleft + idx);
        break;
      }
    }
  }
  gt_suffixsortspace_export_done(sssp);
}
