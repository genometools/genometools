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
#include "core/error.h"
#include "core/encseq.h"
#include "querymatch.h"

typedef int (*GtProcessquerymatch)(void *,
                                   const GtEncseq *,
                                   const GtQuerymatch *,
                                   const GtUchar *query,
                                   unsigned long query_totallength,
                                   GtError *);

typedef void (*GtProcessquerybeforematching)(void *,const char *,
                                             const GtUchar *,unsigned long,
                                             bool);

typedef struct GtMMsearchiterator GtMMsearchiterator;

GtMMsearchiterator *gt_mmsearchiterator_new_complete_olain(
                                    const GtEncseq *dbencseq,
                                    const void *voidsuftab, /* XXX */
                                    unsigned long leftbound,
                                    unsigned long rightbound,
                                    unsigned long itvoffset,
                                    GtReadmode readmode,
                                    const GtUchar *pattern,
                                    unsigned long patternlength);

bool gt_mmsearchiterator_next(unsigned long *dbstart,GtMMsearchiterator *mmsi);

bool gt_mmsearchiterator_isempty(const GtMMsearchiterator *mmsi);

bool gt_mmsearchiterator_identical(const GtMMsearchiterator *mmsi1,
                                   const GtMMsearchiterator *mmsi2);

void gt_mmsearchiterator_delete(GtMMsearchiterator *mmsi);

unsigned long gt_mmsearchiterator_count(const GtMMsearchiterator *mmsi);

int gt_callenumquerymatches(const char *indexname,
                            const GtStrArray *queryfiles,
                            bool findmums,
                            bool forwardstrand,
                            bool reversestrand,
                            unsigned int userdefinedleastlength,
                            GtProcessquerybeforematching
                              processquerybeforematching,
                            GtProcessquerymatch processquerymatch,
                            void *processquerymatchinfo,
                            GtLogger *logger,
                            GtError *err);

int gt_callenumselfmatches(const char *indexname,
                           GtReadmode queryreadmode,
                           unsigned int userdefinedleastlength,
                           GtProcessquerymatch processquerymatch,
                           void *processquerymatchinfo,
                           GtLogger *logger,
                           GtError *err);

int gt_sarrquerysubstringmatch(const GtUchar *dbseq,
                               unsigned long dblen,
                               const GtUchar *query,
                               unsigned long querylen,
                               unsigned int minlength,
                               GtAlphabet *alpha,
                               GtProcessquerymatch processquerymatch,
                               void *processquerymatchinfo,
                               GtLogger *logger,
                               GtError *err);

#endif
