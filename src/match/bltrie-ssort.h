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

#ifndef BLTRIE_SSORT_H
#define BLTRIE_SSORT_H

#include "core/seqpos.h"
#include "encodedsequence.h"
#include "core/readmode.h"

typedef struct Blindtrie Blindtrie;

typedef enum
{
  Ascending,
  Descending,
  Noorder
} Ordertype;

Blindtrie *blindtrie_new(unsigned long numofsuffixes,
                         const GtEncodedsequence *encseq,
                         bool cmpcharbychar,
                         GtEncodedsequenceScanstate *esr1,
                         GtEncodedsequenceScanstate *esr2,
                         GtReadmode readmode);

Seqpos blindtrie_suffixsort(Blindtrie *blindtrie,
                            Seqpos *suffixtable,
                            Seqpos *lcpsubtab,
                            unsigned long numberofsuffixes,
                            Seqpos offset,
                            Seqpos maxdepth,
                            Ordertype ordertype,
                            void *voiddcov,
                            void (*dc_processunsortedrange)(void *,Seqpos *,
                                                            Seqpos *,Seqpos));

void blindtrie_delete(Blindtrie **blindtrie);

#endif
