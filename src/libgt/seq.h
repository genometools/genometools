/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SEQ_H
#define SEQ_H

#include "alpha.h"

typedef struct Seq Seq;

Seq*          seq_new(const char *seq, unsigned long seqlen, Alpha *seqalpha,
                      Env*);
void          seq_set_description(Seq*, const char *desc);
const char*   seq_get_orig(const Seq*); /* not '\0' terminated */
const char*   seq_get_encoded(Seq*);
const Alpha*  seq_get_alpha(const Seq*);
unsigned long seq_length(const Seq*);
void          seq_delete(Seq*, Env*);

#endif
