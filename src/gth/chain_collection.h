/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef CHAIN_COLLECTION_H
#define CHAIN_COLLECTION_H

typedef struct GthChainCollection GthChainCollection;

#include "gth/gthchain.h"

GthChainCollection* gth_chain_collection_new(void);
void                gth_chain_collection_delete(GthChainCollection*);
/* Takes ownership of <chain>. */
void                gth_chain_collection_add(GthChainCollection*,
                                             GthChain *chain);
void                gth_chain_collection_sort(GthChainCollection*);
unsigned long       gth_chain_collection_size(const GthChainCollection*);
GthChain*           gth_chain_collection_get(const GthChainCollection*,
                                             unsigned long);

#endif
