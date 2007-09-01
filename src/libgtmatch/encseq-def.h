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

#ifndef ENCSEQDEF_H
#define ENCSEQDEF_H
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"
#include "symboldef.h"
#include "seqpos-def.h"
#include "alphadef.h"
#include "readmode-def.h"

#define ENCSEQFILESUFFIX     ".esq"

typedef struct Encodedsequence Encodedsequence;
typedef struct Encodedsequencescanstate Encodedsequencescanstate;

typedef struct
{
  Seqpos leftpos,
         rightpos;
} Sequencerange;          /* \Typedef{Sequencerange} */

Seqpos getencseqtotallength(const Encodedsequence *encseq);

Uchar getencodedchar(const Encodedsequence *encseq,Seqpos pos,
                     Readmode readmode);

int flushencseqfile(const Str *indexname,Encodedsequence *encseq,Env *env);

void freeEncodedsequence(Encodedsequence **encseqptr,Env *env);

/*@null@*/ Encodedsequencescanstate *initEncodedsequencescanstate(
                                         const Encodedsequence *encseq,
                                         Readmode readmode,
                                         Env *env);

void freeEncodedsequencescanstate(Encodedsequencescanstate **esr,Env *env);

Uchar sequentialgetencodedchar(const Encodedsequence *encseq,
                               Encodedsequencescanstate *esr,
                               Seqpos pos);

/*@null@*/ Encodedsequence *initencodedseq(bool withrange,
                                           const StrArray *filenametab,
                                           bool plainformat,
                                           const Str *indexname,
                                           Seqpos totallength,
                                           const Specialcharinfo
                                                 *specialcharinfo,
                                           const Alphabet *alphabet,
                                           const char *str_sat,
                                           Env *env);

bool fastspecialranges(const Encodedsequence *encseq);

int overallspecialrangesfast(
                const Encodedsequence *encseq,
                bool moveforward,
                Readmode readmode,
                int(*processrange)(void *,const Encodedsequence *,
                                   const Sequencerange *,Env *),
                void *processinfo,
                Env *env);

int overallspecialranges(const Encodedsequence *encseq,
                         Readmode readmode,
                         int(*processrange)(void *,const Encodedsequence *,
                                            const Sequencerange *,Env *),
                         void *processinfo,
                         Env *env);

#endif
