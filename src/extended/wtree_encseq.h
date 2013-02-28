/*
  Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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
#ifndef WTREE_ENCSEQ_H
#define WTREE_ENCSEQ_H

#include "extended/wtree.h"

#include "core/encseq_api.h"

/* Class <GtWtreeEncseq>, Implements the <GtWtree> interface.
   This implementation represents the sequence part of an encoded sequence */
typedef struct GtWtreeEncseq GtWtreeEncseq;

/* Return a new <GtWtree> object, the implementation beeing of type
   <GtWtreeEncseq>.
   TODO: add documentation */
GtWtree* gt_wtree_encseq_new(GtEncseq *encseq);

/* Mapps a <GtWtreeSymbol> to a <char> according to the encseq it was build with
   */
char gt_wtree_encseq_unmap(GtWtree *wtree, GtWtreeSymbol symbol);
#endif
