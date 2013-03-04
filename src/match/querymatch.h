/*
  Copyright (c) 2007-2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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

#ifndef QUERYMATCH_H
#define QUERYMATCH_H

#include <inttypes.h>
#include "core/error_api.h"
#include "core/readmode.h"
#include "core/encseq.h"

typedef struct GtQuerymatch GtQuerymatch;

GtQuerymatch *gt_querymatch_new(void);

void gt_querymatch_fill(GtQuerymatch *querymatch,
                        unsigned long dblen,
                        unsigned long dbstart,
                        GtReadmode readmode,
                        bool query_as_reversecopy,
                        long score,
                        unsigned long edist,
                        bool selfmatch,
                        uint64_t queryseqnum,
                        unsigned long querylen,
                        unsigned long querystart);

void gt_querymatch_delete(GtQuerymatch *querymatch);

int gt_querymatch_output(void *info,
                         const GtEncseq *encseq,
                         const GtQuerymatch *querymatch,
                         const GtUchar *query,
                         unsigned long query_totallength,
                         GtError *err);

unsigned long gt_querymatch_querylen(const GtQuerymatch *querymatch);

unsigned long gt_querymatch_dbstart(const GtQuerymatch *querymatch);

unsigned long gt_querymatch_querystart(const GtQuerymatch *querymatch);

uint64_t gt_querymatch_queryseqnum(const GtQuerymatch *querymatch);

const GtUchar *gt_querymatch_querysequence(const GtQuerymatch *querymatch);

unsigned long gt_querymatch_querytotallength(const GtQuerymatch *querymatch);

unsigned long gt_querymatch_dbseqnum(const GtEncseq *encseq,
                                     const GtQuerymatch *querymatch);

bool gt_querymatch_queryreverse(const GtQuerymatch *querymatch);

#endif
