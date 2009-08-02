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

#ifndef ALPHABET_H
#define ALPHABET_H

#include <limits.h>
#include "core/str.h"
#include "core/str_array.h"
#include "core/symboldef.h"
#include "core/error_api.h"

#define MAXALPHABETCHARACTER UCHAR_MAX
#define COMPAREOFFSET        (MAXALPHABETCHARACTER + 1)

/*
  the size of the DNA alphabet
*/

#define DNAALPHASIZE        4U

/*
  The following type is for storing alphabets.
*/

typedef struct GtAlphabet GtAlphabet;

/*@null@*/
GtAlphabet*    gt_alphabet_new(bool isdna, bool isprotein,
                               const GtStr *smapfile,
                               const GtStrArray *filenametab, GtError *err);
GtAlphabet*    gt_alphabet_clone(const GtAlphabet *alphabet);
void           gt_alphabet_delete(GtAlphabet *alphabet);
const GtUchar* gt_alphabet_symbolmap(const GtAlphabet *alphabet);
unsigned int   gt_alphabet_num_of_chars(const GtAlphabet *alphabet);
const GtUchar* gt_alphabet_characters(const GtAlphabet *alphabet);
GtUchar        gt_alphabet_wildcard_show(const GtAlphabet *alphabet);
unsigned int   gt_alphabet_bits_per_symbol(const GtAlphabet *alphabet);

void outputalphabet(FILE *fpout,const GtAlphabet *alpha);

void fprintfsymbolstring(FILE *fpout,const GtAlphabet *alpha,
                         const GtUchar *w,unsigned long wlen);

void printfsymbolstring(const GtAlphabet *alpha,const GtUchar *w,
                        unsigned long wlen);

void sprintfsymbolstring(char *buffer,const GtAlphabet *alpha,
                         const GtUchar *w,unsigned long wlen);

void echoprettysymbol(FILE *fpout,const GtAlphabet *alpha,GtUchar currentchar);

GtUchar getprettysymbol(const GtAlphabet *alpha,unsigned int currentchar);

bool isproteinalphabet(const GtAlphabet *alpha);

bool isdnaalphabet(const GtAlphabet *alpha);

#endif
