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

#define REVERSEPOS(TOT,POS) ((TOT) - 1 - (POS))

typedef struct
{
  Seqpos leftpos,
         rightpos;
} Sequencerange;          /* \Typedef{Sequencerange} */

#ifdef INLINEDENCSEQ

typedef struct
{
  Uchar *plainseq;
  Seqpos *totallength;
} Encodedsequence;

#define getencseqtotallength(ENCSEQ) ((ENCSEQ)->totallength)
#define getencodedchar(ENCSEQ,POS,READMODE)\
        (ENCSEQ)->plainseq[POS]

int flushencseqfile(const Str *indexname,Encodedsequence *encseq,Env *env);

#else

typedef struct Encodedsequence Encodedsequence;
typedef struct Encodedsequencescanstate Encodedsequencescanstate;
typedef struct Specialrangeiterator Specialrangeiterator;

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

/*@null@*/ Encodedsequence *files2encodedsequence(bool withrange,
                                                  const StrArray *filenametab,
                                                  bool plainformat,
                                                  Seqpos totallength,
                                                  const Specialcharinfo
                                                        *specialcharinfo,
                                                  const Alphabet *alphabet,
                                                  const char *str_sat,
                                                  Env *env);

/*@null@*/ Encodedsequence *mapencodedsequence(bool withrange,
                                               const Str *indexname,
                                               Seqpos totallength,
                                               const Specialcharinfo
                                                     *specialcharinfo,
                                               const Alphabet *alphabet,
                                               const char *str_sat,
                                               Env *env);

Encodedsequence *plain2encodedsequence(bool withrange,
                                       Specialcharinfo *specialcharinfo,
                                       const Uchar *seq1,
                                       Seqpos len1,
                                       const Uchar *seq2,
                                       unsigned long len2,
                                       const Alphabet *alphabet,
                                       Env *env);

bool hasspecialranges(const Encodedsequence *encseq);

Specialrangeiterator *newspecialrangeiterator(const Encodedsequence *encseq,
                                              bool moveforward,
                                              Env *env);

bool nextspecialrangeiterator(Sequencerange *range,Specialrangeiterator *sri);

void freespecialrangeiterator(Specialrangeiterator **sri,Env *env);

/*@null@*/ char *encseqaccessname(const Encodedsequence *encseq);

#endif

#endif
