/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ALPHA_H
#define ALPHA_H

#include <stdbool.h>
#include <stdio.h>
#include <libgt/env.h>

/* the alphabet class */
typedef struct Alpha Alpha; /* XXX: Alpha -> Alphabet */

Alpha*       alpha_new(Env*);
Alpha*       alpha_new_dna(Env*);
Alpha*       alpha_new_protein(Env*);
Alpha*       alpha_guess(const char *seq, unsigned long seqlen, Env*);
Alpha*       alpha_ref(Alpha*);
/* add the mapping of all given characters to the alphabet, the first character
   is the result of subsequent alpha_decode() calls  */
void         alpha_add_mapping(Alpha*, const char *characters);
char         alpha_decode(const Alpha*, unsigned int);
unsigned int alpha_encode(const Alpha*, char);
void         alpha_decode_seq(const Alpha*, char *out, char *in,
                              unsigned long length); /* in can be == out */
void         alpha_encode_seq(const Alpha*, char *out, char *in,
                              unsigned long length); /* in can be == out */
bool         alpha_is_compatible_with_alpha(const Alpha*, const Alpha*);
unsigned int alpha_size(const Alpha*);
void         alpha_delete(Alpha*, Env*);

#endif
