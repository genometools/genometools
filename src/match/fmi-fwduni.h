/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>

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

#ifndef FMI_FWDUNI_H
#define FMI_FWDUNI_H
#include "core/types_api.h"
#include "core/unused_api.h"

GtUword gt_skfmuniqueforward (const void *genericindex,
                              GT_UNUSED GtUword offset,
                              GT_UNUSED GtUword left,
                              GT_UNUSED GtUword right,
                              GT_UNUSED GtUword *witnessposition,
                              const GtUchar *qstart,
                              const GtUchar *qend);

GtUword gt_skfmmstats (const void *genericindex,
                       GT_UNUSED GtUword offset,
                       GT_UNUSED GtUword left,
                       GT_UNUSED GtUword right,
                       GtUword *witnessposition,
                       const GtUchar *qstart,
                       const GtUchar *qend);

#endif
