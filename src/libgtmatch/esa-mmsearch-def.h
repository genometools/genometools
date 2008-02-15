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

#ifndef ESA_MMSEARCH_DEF_H
#define ESA_MMSEARCH_DEF_H
#include "libgtcore/error.h"
#include "encseq-def.h"

typedef struct MMsearchiterator MMsearchiterator;

MMsearchiterator *newmmsearchiterator(const Encodedsequence *dbencseq,
                                      const Seqpos *suftab,
                                      Seqpos leftbound,
                                      Seqpos rightbound,
                                      Seqpos offset,
                                      Readmode readmode,
                                      const Uchar *pattern,
                                      unsigned long patternlen);

bool nextmmsearchiterator(Seqpos *dbstart,MMsearchiterator *mmsi);

bool isemptymmsearchiterator(const MMsearchiterator *mmsi);

bool identicalmmsearchiterators(const MMsearchiterator *mmsi1,
                                const MMsearchiterator *mmsi2);

void freemmsearchiterator(MMsearchiterator **mmsi);

int runquerysubstringmatch(const Encodedsequence *dbencseq,
                           const Seqpos *suftabpart,
                           Readmode readmode,
                           Seqpos numberofsuffixes,
                           uint64_t unitnum,
                           const Uchar *query,
                           unsigned long querylen,
                           unsigned int minlength,
                           int (*processmaxmatch)(void *,unsigned long,
                                                  Seqpos,uint64_t,
                                                  unsigned long,Error *),
                           void *processmaxmatchinfo,
                           Error *err);

int callenumquerymatches(const Str *indexname,
                         const StrArray *queryfiles,
                         bool echoquery,
                         unsigned int userdefinedleastlength,
                         int (*processmaxmatch)(void *,unsigned long,Seqpos,
                                                uint64_t,unsigned long,Error *),
                         void *processmaxmatchinfo,
                         Verboseinfo *verboseinfo,
                         Error *err);

int sarrquerysubstringmatch(const Uchar *dbseq,
                            Seqpos dblen,
                            const Uchar *query,
                            unsigned long querylen,
                            unsigned int minlength,
                            const Alphabet *alpha,
                            int (*processmaxmatch)(void *,unsigned long,Seqpos,
                                                   uint64_t,unsigned long,
                                                   Error *),
                            void *processmaxmatchinfo,
                            Verboseinfo *verboseinfo,
                            Error *err);

Seqpos countmmsearchiterator(const MMsearchiterator *mmsi);

#endif
