/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SEQ_H
#define SEQ_H

#include <libgtcore/alpha.h>

typedef struct Seq Seq;

/* stores <seq> pointer */
Seq*          seq_new(const char *seq, unsigned long seqlen, Alpha *seqalpha,
                      Env*);
/* takes ownership of <seq> */
Seq*          seq_new_own(char *seq, unsigned long seqlen,
                          Alpha *seqalpha, Env*);
/* stores <desc> pointer */
void          seq_set_description(Seq*, const char *desc);
/* takes ownership of <desc> */
void          seq_set_description_own(Seq*, char *desc, Env*);
void          seq_set_description(Seq*, const char *desc);
const char*   seq_get_description(Seq*);
const char*   seq_get_orig(const Seq*); /* not '\0' terminated */
const char*   seq_get_encoded(Seq*, Env*);
const Alpha*  seq_get_alpha(const Seq*);
unsigned long seq_length(const Seq*);
void          seq_delete(Seq*, Env*);

#endif
