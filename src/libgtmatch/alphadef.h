/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ALPHADEF_H
#define ALPHADEF_H

#include "libgtcore/str.h"
#include "libgtcore/strarray.h"
#include "symboldef.h"

/*
  The following type is for storing alphabets.
*/

typedef struct Alphabet Alphabet;

/*@null@*/ Alphabet *assigninputalphabet(bool isdna,
                                         bool isprotein,
                                         const Str *smapfile,
                                         const StrArray *filenametab,
                                         Env *env);

const Uchar *getsymbolmapAlphabet(const Alphabet *alpha);

unsigned int getnumofcharsAlphabet(const Alphabet *alpha);

unsigned int getmapsizeAlphabet(const Alphabet *alpha);

const Uchar *getcharactersAlphabet(const Alphabet *alpha);

Uchar *copycharactersAlphabet(const Alphabet *alpha,Env *env);

void freeAlphabet(Alphabet **alpha,Env *env);

void outputalphabet(FILE *fpout,const Alphabet *alpha);

void showsymbolstringgeneric(FILE *fpout,const Alphabet *alpha,
                             const Uchar *w,unsigned long wlen);

void showsymbolstring(const Alphabet *alpha,const Uchar *w,unsigned long wlen);

bool isproteinalphabet(const Alphabet *alpha);

bool isdnaalphabet(const Alphabet *alpha);

#endif
