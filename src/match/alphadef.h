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

#include <limits.h>
#include "core/str.h"
#include "core/str_array.h"
#include "core/symboldef.h"

#define MAXALPHABETCHARACTER UCHAR_MAX

/*
  the size of the DNA alphabet
*/

#define DNAALPHASIZE        4U
#define PROTEINALPHASIZE   20U
#define BITSFORAMINOACID    5U

/*
  The following type is for storing alphabets.
*/

typedef struct SfxAlphabet SfxAlphabet;

/*@null@*/ SfxAlphabet *assigninputalphabet(bool isdna,
                                            bool isprotein,
                                            const GtStr *smapfile,
                                            const GtStrArray *filenametab,
                                            GtError *err);

SfxAlphabet *gt_copyAlphabet(const SfxAlphabet *alpha2);

const Uchar *getsymbolmapAlphabet(const SfxAlphabet *alpha);

unsigned int getnumofcharsAlphabet(const SfxAlphabet *alpha);

const Uchar *getcharactersAlphabet(const SfxAlphabet *alpha);

Uchar getwildcardshowAlphabet(const SfxAlphabet *alpha);

void freeSfxAlphabet(SfxAlphabet **alpha);

void outputalphabet(FILE *fpout,const SfxAlphabet *alpha);

void fprintfsymbolstring(FILE *fpout,const SfxAlphabet *alpha,
                         const Uchar *w,unsigned long wlen);

void printfsymbolstring(const SfxAlphabet *alpha,const Uchar *w,
                        unsigned long wlen);

void sprintfsymbolstring(char *buffer,const SfxAlphabet *alpha,
                         const Uchar *w,unsigned long wlen);

void echoprettysymbol(FILE *fpout,const SfxAlphabet *alpha,Uchar currentchar);

Uchar getprettysymbol(const SfxAlphabet *alpha,unsigned int currentchar);

bool isproteinalphabet(const SfxAlphabet *alpha);

bool isdnaalphabet(const SfxAlphabet *alpha);

#endif
