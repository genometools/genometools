/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SPLITTER_H
#define SPLITTER_H

typedef struct Splitter Splitter;

Splitter*     splitter_new(void);
void          splitter_split(Splitter*, char *, unsigned long, char delimiter);
char**        splitter_get_tokens(Splitter*);
char*         splitter_get_token(Splitter*, unsigned long);
void          splitter_reset(Splitter*);
unsigned long splitter_size(Splitter*);
int           splitter_unit_test(void);
void          splitter_free(Splitter*);

#endif
