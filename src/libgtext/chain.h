/*
  Copyright (c)       2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004       Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2004, 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef CHAIN_H
#define CHAIN_H

#include "libgtcore/error.h"

typedef struct Chain Chain;

Chain*        chain_new(void);
void          chain_reset(Chain*);
long          chain_get_score(const Chain*);
void          chain_set_score(Chain*, long);
void          chain_add_fragnum(Chain*, unsigned long fragnum);
void          chain_set_fragnum(Chain*, unsigned long idx,
                                        unsigned long fragnum);
unsigned long chain_get_fragnum(const Chain*, unsigned long idx);
unsigned long chain_size(const Chain*);
void          chain_delete(Chain*);

#endif
