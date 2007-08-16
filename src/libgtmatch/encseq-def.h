/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
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

int overallspecialrangesforward(
                const Encodedsequence *encseq,
                Readmode readmode,
                int(*process)(void *,const Sequencerange *,Env *),
                void *processinfo,
                Env *env);

int overallspecialrangesbackward2(
                const Encodedsequence *encseq,
                Readmode readmode,
                int(*process)(void *,const Sequencerange *,Env *),
                void *processinfo,
                Env *env);

int overallspecialranges(const Encodedsequence *encseq,
                         Readmode readmode,
                         int(*process)(void *,const Sequencerange *,Env *),
                         void *processinfo,
                         Env *env);

#endif
