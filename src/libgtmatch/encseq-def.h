/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ENCSEQDEF_H
#define ENCSEQDEF_H
#include "types.h"
#include "alphadef.h"

typedef struct _Encodedsequence Encodedsequence;
typedef struct _Encodedsequencescanstate Encodedsequencescanstate;

Uint64 getencseqtotallength(const Encodedsequence *encseq);

Uchar getencodedchar(const Encodedsequence *encseq,Uint pos);

Uchar getencodedchar64(const Encodedsequence *encseq,Uint64 pos);

int flushencseqfile(const char *indexname,Encodedsequence *encseq,Env *env);

void freeEncodedsequence(Encodedsequence **encseqptr,Env *env);

/*@null@*/ Encodedsequencescanstate *initEncodedsequencescanstate(
                                         const Encodedsequence *encseq,
                                         Env *env);

void freeEncodedsequencescanstate(Encodedsequencescanstate **esr,Env *env);

Uchar sequentialgetencodedchar64(const Encodedsequence *encseq,
                                 Encodedsequencescanstate *esr,
                                 Uint64 pos);

int overallspecialranges(const Encodedsequence *encseq,
                         int(*process)(void *,const PairUint64 *,Env *),
                         void *processinfo,
                         Env *env);

/*@null@*/ Encodedsequence *initencodedseq(bool withrange,
                                           const char **filenametab,
                                           unsigned int numoffiles,
                                           const char *indexname,
                                           Uint64 totallength,
                                           const Specialcharinfo
                                                 *specialcharinfo,
                                           const Alphabet *alphabet,
                                           const char *str_sat,
                                           Env *env);

#endif
