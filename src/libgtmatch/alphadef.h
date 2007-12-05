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

#ifndef ALPHADEF_H
#define ALPHADEF_H

#include "libgtcore/str.h"
#include "libgtcore/strarray.h"
#include "libgtcore/symboldef.h"

/*
  The following type is for storing alphabets.
*/

typedef struct Alphabet Alphabet;

/*@null@*/ Alphabet *assigninputalphabet(bool isdna,
                                         bool isprotein,
                                         const Str *smapfile,
                                         const StrArray *filenametab,
                                         Error *err);

const Uchar *getsymbolmapAlphabet(const Alphabet *alpha);

unsigned int getnumofcharsAlphabet(const Alphabet *alpha);

unsigned int getmapsizeAlphabet(const Alphabet *alpha);

const Uchar *getcharactersAlphabet(const Alphabet *alpha);

void freeAlphabet(Alphabet **alpha);

void outputalphabet(FILE *fpout,const Alphabet *alpha);

void showsymbolstringgeneric(FILE *fpout,const Alphabet *alpha,
                             const Uchar *w,unsigned long wlen);

void showsymbolstring(const Alphabet *alpha,const Uchar *w,unsigned long wlen);

void echoprettysymbol(FILE *fpout,const Alphabet *alpha,Uchar currentchar);

Uchar getprettysymbol(const Alphabet *alpha,unsigned int currentchar);

bool isproteinalphabet(const Alphabet *alpha);

bool isdnaalphabet(const Alphabet *alpha);

#endif
