/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ENCSEQDEF_H
#define ENCSEQDEF_H
#include "types.h"
#include "alphadef.h"
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"

typedef struct _Encodedsequence Encodedsequence;
typedef struct _Encodedsequencescanstate Encodedsequencescanstate;

typedef struct
{
  Seqpos leftpos,
         rightpos;
} Sequencerange;          /* \Typedef{Sequencerange} */

Seqpos getencseqtotallength(const Encodedsequence *encseq);

Uchar getencodedchar(const Encodedsequence *encseq,Seqpos pos);

int flushencseqfile(const Str *indexname,Encodedsequence *encseq,Env *env);

void freeEncodedsequence(Encodedsequence **encseqptr,Env *env);

/*@null@*/ Encodedsequencescanstate *initEncodedsequencescanstate(
                                         const Encodedsequence *encseq,
                                         Env *env);

void freeEncodedsequencescanstate(Encodedsequencescanstate **esr,Env *env);

Uchar sequentialgetencodedchar(const Encodedsequence *encseq,
                               Encodedsequencescanstate *esr,
                               Seqpos pos);

int overallspecialranges(const Encodedsequence *encseq,
                         int(*process)(void *,const Sequencerange *,Env *),
                         void *processinfo,
                         Env *env);

/*@null@*/ Encodedsequence *initencodedseq(bool withrange,
                                           const StrArray *filenametab,
                                           const Str *indexname,
                                           Seqpos totallength,
                                           const Specialcharinfo
                                                 *specialcharinfo,
                                           const Alphabet *alphabet,
                                           bool plainformat,
                                           const char *str_sat,
                                           Env *env);

#endif
