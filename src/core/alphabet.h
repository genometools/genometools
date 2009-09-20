/*
  Copyright (c) 2007-2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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

#define GT_MAXALPHABETCHARACTER UCHAR_MAX
#define GT_COMPAREOFFSET        (GT_MAXALPHABETCHARACTER + 1)

/*
  the size of the DNA alphabet
*/

#define GT_DNAALPHASIZE        4U

/*
  The following type is for storing alphabets.
*/

typedef struct GtAlphabet GtAlphabet;

/*@null@*/
GtAlphabet*    gt_alphabet_new(bool isdna, bool isprotein,
                               const GtStr *smapfile,
                               const GtStrArray *filenametab, GtError *err);
/* Return a <GtAlphabet*> object which represents a DNA alphabet. */
GtAlphabet*    gt_alphabet_new_dna(void);
/* Return a <GtAlphabet*> object which represents a protein alphabet. */
GtAlphabet*    gt_alphabet_new_protein(void);
/* Return an empty <GtAlphabet*> object. */
GtAlphabet*    gt_alphabet_new_empty(void);
/* Try to guess which type the given <sequence> with <length> has (DNA or
   protein) and return an according <GtAlphabet*> object. */
GtAlphabet*    gt_alphabet_guess(const char *sequence, unsigned long seqlen);
/* Return a clone of <alphabet>. */
GtAlphabet*    gt_alphabet_clone(const GtAlphabet *alphabet);
/* Increase the reference count for <alphabet> and return it. */
GtAlphabet*    gt_alphabet_ref(GtAlphabet *alphabet);
/* Decrease the reference count for <alphabet> or delete it, if this was the
   last reference. */
void           gt_alphabet_delete(GtAlphabet *alphabet);
/* Add the mapping of all given <characters> to the given <alphabet>. The first
   character is the result of subsequent <gt_alphabet_decode()> calls. */
void           gt_alphabet_add_mapping(GtAlphabet *alphabet,
                                       const char *characters);
/* Add <wildcard> to the <alphabet>. */
void           gt_alphabet_add_wildcard(GtAlphabet *alphabet, char wildcard);
const GtUchar* gt_alphabet_symbolmap(const GtAlphabet *alphabet);
unsigned int   gt_alphabet_num_of_chars(const GtAlphabet *alphabet);
unsigned int   gt_alphabet_size(const GtAlphabet *alphabet);
const GtUchar* gt_alphabet_characters(const GtAlphabet *alphabet);
GtUchar        gt_alphabet_wildcard_show(const GtAlphabet *alphabet);
unsigned int   gt_alphabet_bits_per_symbol(const GtAlphabet *alphabet);
void           gt_alphabet_output(const GtAlphabet *alphabet, FILE *fpout);
/* Suppose the string <w> of length <wlen> was transformed according to the
   <alphabet>. The following method shows each character in <w> as the
   characters specified in the transformation. The output is written to the
   given file pointer <fpout>. */
void           gt_alphabet_fprintf_symbolstring(const GtAlphabet *alphabet,
                                                FILE *fpout, const GtUchar *w,
                                                unsigned long wlen);
/* Analog to <gt_alphabet_fprintf_symbolstring()> showing the output on
   <stdout>. */
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
/* The following method checks if the given <alphabet> is the protein
   alphabet with the aminoacids A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S,
   T, V, W, Y written in lower or upper case and returns <true>, if this is the
   case (<false> otherwise). */
bool           gt_alphabet_is_protein(const GtAlphabet *alphabet);
/* The following method checks if the given alphabet is the DNA alphabet with
   the bases A, C, G, T written in lower or upper case and returns <true>, if
   this is the case (<false> otherwise). */
bool           gt_alphabet_is_dna(const GtAlphabet *alphabet);
bool           gt_alphabet_valid_input(const GtAlphabet *alphabet, char c);
/* Encode character <c> with given <alphabet>.
   <c> has to be encodable with the given <alphabet>! */
GtUchar        gt_alphabet_encode(const GtAlphabet *alphabet, char c);
/* Decode character <c> with given <alphabet>. */
char           gt_alphabet_decode(const GtAlphabet *alphabet, GtUchar c);
/* Encode sequence <in> of given <length> with <alphabet> and store the result
   in <out>. <in> has to be encodable with the given <alphabet>! */
void           gt_alphabet_encode_seq(const GtAlphabet *alphabet, GtUchar *out,
                                      const char *in, unsigned long length);

#endif
