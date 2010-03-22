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

#include "splititv.h"
#include "intcode-def.h"

typedef struct Pckbuckettable Pckbuckettable;

void pckbuckettable_free(Pckbuckettable *pckbt);

Pckbuckettable *pckbuckettable_new(const void *voidbwtseq,
                                   unsigned int numofchars,
                                   unsigned long totallength,
                                   unsigned int maxdepth);

int pckbucket2file(const GtStr *indexname,const Pckbuckettable *pckbuckettable,
                   GtError *err);

bool pckbuckettableexists(const GtStr *indexname);

Pckbuckettable *mappckbuckettable(const GtStr *indexname,
                                  unsigned int numofchars,
                                  GtError *err);

void enumlowlevelchildintervals(GtArrayBoundswithchar *bwci,
                                const Pckbuckettable *pcktb,
                                GtCodetype parentcode,
                                unsigned long childdepth);

unsigned int pcktb2maxdepth(const Pckbuckettable *pckbuckettable);

const void *pcktb2mbtab(const Pckbuckettable *pckbuckettable);

#endif
