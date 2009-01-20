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

#include "core/error.h"

typedef struct GtChain GtChain;

GtChain*      gt_chain_new(void);
void          gt_chain_reset(GtChain*);
long          gt_chain_get_score(const GtChain*);
void          gt_chain_set_score(GtChain*, long);
void          gt_chain_add_fragnum(GtChain*, unsigned long fragnum);
void          gt_chain_set_fragnum(GtChain*, unsigned long idx,
                                             unsigned long fragnum);
unsigned long gt_chain_get_fragnum(const GtChain*, unsigned long idx);
unsigned long gt_chain_size(const GtChain*);
void          gt_chain_delete(GtChain*);

#endif
