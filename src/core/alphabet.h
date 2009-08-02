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
void           gt_alphabet_output(const GtAlphabet *alphabet, FILE *fpout);
void           gt_alphabet_fprintf_symbolstring(const GtAlphabet *alphabet,
                                                FILE *fpout, const GtUchar *w,
                                                unsigned long wlen);
void           gt_alphabet_printf_symbolstring(const GtAlphabet *alphabet,
                                               const GtUchar *w,
                                               unsigned long wlen);
void           gt_alphabet_sprintf_symbolstring(const GtAlphabet *alphabet,
                                                char *buffer, const GtUchar *w,
                                                unsigned long wlen);
void           gt_alphabet_echo_pretty_symbol(const GtAlphabet *alphabet,
                                              FILE *fpout, GtUchar currentchar);
GtUchar        gt_alphabet_pretty_symbol(const GtAlphabet *alphabet,
                                         unsigned int currentchar);
bool           gt_alphabet_is_protein(const GtAlphabet *alphabet);
bool           gt_alphabet_is_dna(const GtAlphabet *alphabet);

#endif
