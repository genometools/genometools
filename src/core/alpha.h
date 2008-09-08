/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef ALPHA_H
#define ALPHA_H

#include <stdbool.h>
#include <stdio.h>

/* the alphabet class */
typedef struct GT_Alpha GT_Alpha; /* XXX: GT_Alpha -> GT_Alphabet */

GT_Alpha*    gt_alpha_new(void);
GT_Alpha*    gt_alpha_new_dna(void);
GT_Alpha*    gt_alpha_new_protein(void);
GT_Alpha*    gt_alpha_guess(const char *seq, unsigned long seqlen);
GT_Alpha*    gt_alpha_ref(GT_Alpha*);
/* add the mapping of all given characters to the alphabet, the first
   character is the result of subsequent gt_alpha_decode() calls  */
void         gt_alpha_add_mapping(GT_Alpha*, const char *characters);
char         gt_alpha_decode(const GT_Alpha*, unsigned int);
unsigned int gt_alpha_encode(const GT_Alpha*, char);
void         gt_alpha_decode_seq(const GT_Alpha*, char *out, char *in,
                                 unsigned long length); /* in can be == out */
void         gt_alpha_encode_seq(const GT_Alpha*, char *out, char *in,
                                 unsigned long length); /* in can be == out */
bool         gt_alpha_char_is_valid(const GT_Alpha*, char);
bool         gt_alpha_is_compatible_with_alpha(const GT_Alpha*,
                                               const GT_Alpha*);
unsigned int gt_alpha_size(const GT_Alpha*);
void         gt_alpha_delete(GT_Alpha*);

#endif
