/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SPLITTER_H
#define SPLITTER_H

typedef struct Splitter Splitter;

Splitter*     splitter_new(void);

/* split 'string' of given 'length' into tokens delimited by 'delimiter'.
   'string' is modified in the splitting process! */
void          splitter_split(Splitter*, char *string, unsigned long length,
                             char delimiter);

/* get all tokens */
char**        splitter_get_tokens(Splitter*);

/* get token with number 'token_num' */
char*         splitter_get_token(Splitter*, unsigned long token_num);

/* reset the splitter */
void          splitter_reset(Splitter*);

/* returns the number of tokens */
unsigned long splitter_size(Splitter*);

int           splitter_unit_test(void);
void          splitter_free(Splitter*);

#endif
