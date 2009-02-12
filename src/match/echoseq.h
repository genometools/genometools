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

#ifndef ECHOSEQ_H
#define ECHOSEQ_H
#include <stdio.h>
#include "core/symboldef.h"
#include "core/str_array_api.h"
#include "core/error_api.h"
#include "alphadef.h"
#include "readmode-def.h"
#include "seqpos-def.h"
#include "encseq-def.h"

void symbolstring2fasta(FILE *fpout,
                        const char *desc,
                        const SfxAlphabet *alpha,
                        const Uchar *w,
                        unsigned long wlen,
                        unsigned long width);

void encseq2symbolstring(FILE *fpout,
                         const Encodedsequence *encseq,
                         Readmode readmode,
                         Seqpos start,
                         Seqpos wlen,
                         unsigned long width);

void fprintfencseq(FILE *fpout,
                   const Encodedsequence *encseq,
                   Seqpos start,
                   Seqpos wlen);

void encseq2fastaoutput(FILE *fpout,
                        const char *desc,
                        const Encodedsequence *encseq,
                        Readmode readmode,
                        Seqpos start,
                        Seqpos wlen,
                        unsigned long width);

int echodescriptionandsequence(const GtStrArray *filenametab,GtError *err);

#endif
