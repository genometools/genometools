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
#include "encseq-def.h"

typedef int (*Processquerymatch)(void *,unsigned long,Seqpos,Readmode,
                                 uint64_t,unsigned long,GtError *);

typedef struct MMsearchiterator MMsearchiterator;

MMsearchiterator *newmmsearchiteratorcomplete_plain(
                                    const Encodedsequence *dbencseq,
                                    const Seqpos *suftab,
                                    Seqpos leftbound,
                                    Seqpos rightbound,
                                    Seqpos itvoffset,
                                    Readmode readmode,
                                    const GtUchar *pattern,
                                    unsigned long patternlength);

bool nextmmsearchiterator(Seqpos *dbstart,MMsearchiterator *mmsi);

bool isemptymmsearchiterator(const MMsearchiterator *mmsi);

bool identicalmmsearchiterators(const MMsearchiterator *mmsi1,
                                const MMsearchiterator *mmsi2);

void freemmsearchiterator(MMsearchiterator **mmsi);

int callenumquerymatches(const GtStr *indexname,
                         const GtStrArray *queryfiles,
                         bool echoquery,
                         unsigned int userdefinedleastlength,
                         Processquerymatch processquerymatch,
                         void *processquerymatchinfo,
                         Verboseinfo *verboseinfo,
                         GtError *err);

int callenumselfmatches(const GtStr *indexname,
                        Readmode queryreadmode,
                        unsigned int userdefinedleastlength,
                        Processquerymatch processquerymatch,
                        void *processquerymatchinfo,
                        Verboseinfo *verboseinfo,
                        GtError *err);

int sarrquerysubstringmatch(const GtUchar *dbseq,
                            Seqpos dblen,
                            const GtUchar *query,
                            unsigned long querylen,
                            unsigned int minlength,
                            const GtAlphabet *alpha,
                            Processquerymatch processquerymatch,
                            void *processquerymatchinfo,
                            Verboseinfo *verboseinfo,
                            GtError *err);

Seqpos countmmsearchiterator(const MMsearchiterator *mmsi);

#endif
