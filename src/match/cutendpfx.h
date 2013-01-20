/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef CUTENDPFX_H
#define CUTENDPFX_H

#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include "core/types_api.h"
#include "lcpinterval.h"
#include "bcktab.h"

typedef struct Bucketenumerator Bucketenumerator;

Bucketenumerator *gt_newbucketenumerator(const GtBcktab *bcktab,
                                         unsigned int prefixlength,
                                         const GtUchar *demandprefix,
                                         unsigned int demandprefixlength);

bool gt_nextbucketenumerator(Lcpinterval *itv,
                             Bucketenumerator *bucketenumerator);

void gt_freebucketenumerator(Bucketenumerator *bucketenumerator);

#endif
