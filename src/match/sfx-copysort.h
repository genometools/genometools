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

#ifndef SFX_COPYSORT_H
#define SFX_COPYSORT_H

#include "core/codetype.h"
#include "core/encseq.h"
#include "core/logger.h"
#include "bcktab.h"
#include "sfx-suffixgetset.h"

typedef struct GtBucketspec2 GtBucketspec2;

GtBucketspec2 *gt_copysort_new(const GtBcktab *bcktab,
                               const GtEncseq *encseq,
                               GtReadmode readmode,
                               unsigned long partwidth,
                               unsigned int numofchars);

void gt_copysort_derivesorting(const GtBucketspec2 *bucketspec2,
                               GtSuffixsortspace *suffixsortspace,
                               GtLogger *logger);

bool gt_copysort_checkhardwork(const GtBucketspec2 *bucketspec2,
                               GtCodetype code);

void gt_copysort_delete(GtBucketspec2 *bucketspec2);

#endif
