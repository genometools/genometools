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

typedef int (*Processquerymatch)(void *,
                                 const GtEncseq *,
                                 const Querymatch *,
                                 GtError *);

typedef struct MMsearchiterator MMsearchiterator;

MMsearchiterator *gt_newmmsearchiteratorcomplete_plain(
                                    const GtEncseq *dbencseq,
                                    const unsigned long *suftab,
                                    unsigned long leftbound,
                                    unsigned long rightbound,
                                    unsigned long itvoffset,
                                    GtReadmode readmode,
                                    const GtUchar *pattern,
                                    unsigned long patternlength);

bool gt_nextmmsearchiterator(unsigned long *dbstart,MMsearchiterator *mmsi);

bool gt_isemptymmsearchiterator(const MMsearchiterator *mmsi);

bool gt_identicalmmsearchiterators(const MMsearchiterator *mmsi1,
                                const MMsearchiterator *mmsi2);

void gt_freemmsearchiterator(MMsearchiterator **mmsi);

int gt_callenumquerymatches(const GtStr *indexname,
                         const GtStrArray *queryfiles,
                         bool echoquery,
                         unsigned int userdefinedleastlength,
                         Processquerymatch processquerymatch,
                         void *processquerymatchinfo,
                         GtLogger *logger,
                         GtError *err);

int gt_callenumselfmatches(const GtStr *indexname,
                        GtReadmode queryreadmode,
                        unsigned int userdefinedleastlength,
                        Processquerymatch processquerymatch,
                        void *processquerymatchinfo,
                        GtLogger *logger,
                        GtError *err);

int gt_sarrquerysubstringmatch(const GtUchar *dbseq,
                            unsigned long dblen,
                            const GtUchar *query,
                            unsigned long querylen,
                            unsigned int minlength,
                            GtAlphabet *alpha,
                            Processquerymatch processquerymatch,
                            void *processquerymatchinfo,
                            GtLogger *logger,
                            GtError *err);

unsigned long gt_countmmsearchiterator(const MMsearchiterator *mmsi);

#endif
