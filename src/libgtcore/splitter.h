/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef SPLITTER_H
#define SPLITTER_H

#include "libgtcore/error.h"

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

int           splitter_unit_test(Error*);
void          splitter_delete(Splitter*);

#endif
