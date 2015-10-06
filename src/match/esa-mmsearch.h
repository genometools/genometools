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

#ifndef ESA_MMSEARCH_H
#define ESA_MMSEARCH_H
#include "core/error_api.h"
#include "core/encseq.h"
#include "sarr-def.h"
#include "querymatch.h"

typedef void (*GtProcessquerymatch)(void *,const GtQuerymatch *);

typedef struct GtMMsearchiterator GtMMsearchiterator;

GtMMsearchiterator *gt_mmsearchiterator_new_complete_plain(
                                    const GtEncseq *dbencseq,
                                    const void *voidsuftab, /* XXX */
                                    GtUword leftbound,
                                    GtUword rightbound,
                                    GtUword itvoffset,
                                    GtReadmode readmode,
                                    const GtUchar *pattern,
                                    GtUword patternlength);

bool gt_mmsearchiterator_next(GtUword *dbstart,GtMMsearchiterator *mmsi);

bool gt_mmsearchiterator_isempty(const GtMMsearchiterator *mmsi);

bool gt_mmsearchiterator_identical(const GtMMsearchiterator *mmsi1,
                                   const GtMMsearchiterator *mmsi2);

void gt_mmsearchiterator_delete(GtMMsearchiterator *mmsi);

GtUword gt_mmsearchiterator_count(const GtMMsearchiterator *mmsi);

int gt_sarrquerysubstringmatch(const GtUchar *dbseq,
                               GtUword dblen,
                               const GtUchar *query,
                               GtUword querylen,
                               unsigned int minlength,
                               GtAlphabet *alpha,
                               GtProcessquerymatch processquerymatch,
                               void *processquerymatchinfo,
                               GtLogger *logger,
                               GtError *err);

typedef struct GtQuerysubstringmatchiterator GtQuerysubstringmatchiterator;

GtQuerysubstringmatchiterator *gt_querysubstringmatchiterator_new(
                                     const GtEncseq *dbencseq,
                                     GtUword totallength,
                                     const ESASuffixptr *suftabpart,
                                     GtReadmode readmode,
                                     GtUword numberofsuffixes,
                                     const GtStrArray *query_files,
                                     const GtEncseq *query_encseq,
                                     GtReadmode query_readmode,
                                     unsigned int userdefinedleastlength,
                                     GtError *err);

void gt_querysubstringmatchiterator_delete(GtQuerysubstringmatchiterator *qsmi);

GtUword gt_querysubstringmatchiterator_dbstart(
                      const GtQuerysubstringmatchiterator *qsmi);

GtUword gt_querysubstringmatchiterator_querystart(
                      const GtQuerysubstringmatchiterator *qsmi);

GtUword gt_querysubstringmatchiterator_matchlength(
                      const GtQuerysubstringmatchiterator *qsmi);

uint64_t gt_querysubstringmatchiterator_queryunitnum(
                      const GtQuerysubstringmatchiterator *qsmi);

GtUword gt_querysubstringmatchiterator_query_seqlen(
                      const GtQuerysubstringmatchiterator *qsmi);

const GtUchar *gt_querysubstringmatchiterator_query(
                      const GtQuerysubstringmatchiterator *qsmi);

int gt_querysubstringmatchiterator_next(GtQuerysubstringmatchiterator *qsmi,
                                        GtError *err);

#endif
