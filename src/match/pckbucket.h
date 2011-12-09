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

#ifndef PCKBUCKET_H
#define PCKBUCKET_H

#include "core/str.h"
#include "core/error.h"
#include "core/codetype.h"
#include "eis-voiditf.h"
#include "splititv.h"

typedef struct Pckbuckettable Pckbuckettable;

void gt_pckbuckettable_delete(Pckbuckettable *pckbt);

Pckbuckettable *gt_pckbuckettable_new(const FMindex *fmindex,
                                      unsigned int numofchars,
                                      unsigned long totallength,
                                      unsigned int maxdepth);

int gt_pckbuckettable_2file(const char *indexname,
                            const Pckbuckettable *pckbuckettable,
                            GtError *err);

bool gt_pckbuckettable_exists(const char *indexname);

Pckbuckettable *gt_pckbuckettable_map(const char *indexname,
                                      unsigned int numofchars,
                                      GtError *err);

unsigned int gt_pckbuckettable_maxdepth_get(const Pckbuckettable
                                             *pckbuckettable);

const void *gt_pckbuckettable_mbtab_get(const Pckbuckettable *pckbuckettable);

unsigned int gt_pckbuckettable_numofchars_get(
                       const Pckbuckettable *pckbuckettable);

#endif
