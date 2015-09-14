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
#include "core/unused_api.h"
#include "querymatch-align.h"

typedef struct GtQuerymatch GtQuerymatch;

GtQuerymatch *gt_querymatch_new(GtQuerymatchoutoptions *querymatchoutoptions);

void gt_querymatch_init(GtQuerymatch *querymatch,
                        GtUword dblen,
                        GtUword dbstart,
                        GtUword dbseqnum,
                        GtUword dbstart_relative,
                        GtReadmode readmode,
                        bool query_as_reversecopy,
                        GtWord score,
                        GtUword distance,
                        bool selfmatch,
                        uint64_t queryseqnum,
                        GtUword querylen,
                        GtUword querystart,
                        GtUword query_totallength);

void gt_querymatch_delete(GtQuerymatch *querymatch);

int gt_querymatch_output(GT_UNUSED void *info,
                         GT_UNUSED const GtEncseq *encseq,
                         const GtQuerymatch *querymatch,
                         GT_UNUSED const GtUchar *query,
                         GT_UNUSED GtUword query_totallength,
                         GT_UNUSED GtError *err);

bool gt_querymatch_complete(GtQuerymatch *querymatchptr,
                            GtUword dblen,
                            GtUword dbstart,
                            GtUword dbseqnum,
                            GtUword dbstart_relative,
                            GtReadmode readmode,
                            bool query_as_reversecopy,
                            GtWord score,
                            GtUword distance,
                            bool selfmatch,
                            uint64_t queryseqnum,
                            GtUword querylen,
                            GtUword querystart,
                            const GtEncseq *encseq,
                            const GtUchar *query,
                            GtUword query_totallength,
                            GtUword seedpos1,
                            GtUword seedpos2,
                            GtUword seedlen,
                            bool greedyextension);

GtUword gt_querymatch_querylen(const GtQuerymatch *querymatch);

GtUword gt_querymatch_dbstart(const GtQuerymatch *querymatch);

GtUword gt_querymatch_querystart(const GtQuerymatch *querymatch);

uint64_t gt_querymatch_queryseqnum(const GtQuerymatch *querymatch);

const GtUchar *gt_querymatch_querysequence(const GtQuerymatch *querymatch);

GtUword gt_querymatch_querytotallength(const GtQuerymatch *querymatch);

GtUword gt_querymatch_dbseqnum(const GtQuerymatch *querymatch);

bool gt_querymatch_queryreverse(const GtQuerymatch *querymatch);

double gt_querymatch_error_rate(GtUword distance,GtUword alignedlen);

void gt_querymatch_prettyprint(const GtQuerymatch *querymatch);

GtWord gt_querymatch_distance2score(GtUword distance,GtUword alignedlen);

/* Returns true, iff the given seed start position is below the querymatch. */
bool gt_querymatch_checkoverlap(const GtQuerymatch *querymatch,
                                GtUword start_relative);
#endif
