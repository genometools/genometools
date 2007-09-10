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

#ifndef ESA_SEQREAD_H
#define ESA_SEQREAD_H
#include <stdbool.h>
#include "libgtcore/str.h"
#include "libgtcore/env.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "sarr-def.h"

typedef enum
{
  SEQ_mappedboth,
  SEQ_scan,
  SEQ_suftabfrommemory
} Sequentialaccesstype;

typedef struct Sequentialsuffixarrayreader Sequentialsuffixarrayreader;

Sequentialsuffixarrayreader *newSequentialsuffixarrayreaderfromfile(
                                        const Str *indexname,
                                        unsigned int demand,
                                        Sequentialaccesstype seqactype,
                                        Env *env);

Sequentialsuffixarrayreader *newSequentialsuffixarrayreaderfromRAM(
                                        const Encodedsequence *encseq,
                                        const Seqpos *suftab,
                                        Seqpos numberofsuffixes,
                                        Readmode readmode,
                                        Env *env);

void freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar,
                                     Env *env);

int nextSequentiallcpvalue(Seqpos *currentlcp,
                           Sequentialsuffixarrayreader *ssar,
                           Env *env);

int nextSequentialsuftabvalue(Seqpos *currentsuffix,
                              Sequentialsuffixarrayreader *ssar,
                              Env *env);

const Encodedsequence *encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *sarr);

Readmode readmodeSequentialsuffixarrayreader(
			  const Sequentialsuffixarrayreader *sarr);

const Alphabet *alphabetSequentialsuffixarrayreader(
			  const Sequentialsuffixarrayreader *sarr);

unsigned long numofdbsequencesSequentialsuffixarrayreader(
                    const Sequentialsuffixarrayreader *sarr);

#endif
