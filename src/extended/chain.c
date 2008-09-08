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

#include "core/array.h"
#include "core/ma.h"
#include "extended/chain.h"

struct Chain {
  GT_Array *fragments;
  long score;
};

Chain* chain_new(void)
{
  Chain *chain;
  chain = gt_calloc(1, sizeof *chain);
  chain->fragments = gt_array_new(sizeof (unsigned long));
  return chain;
}

void chain_reset(Chain *chain)
{
  assert(chain);
  chain->score = 0;
  gt_array_reset(chain->fragments);
}

long chain_get_score(const Chain *chain)
{
  assert(chain);
  return chain->score;
}

void chain_set_score(Chain *chain, long score)
{
  assert(chain);
  chain->score = score;
}

void chain_add_fragnum(Chain *chain, unsigned long fragnum)
{
  assert(chain);
  gt_array_add(chain->fragments, fragnum);
}

void chain_set_fragnum(Chain *chain, unsigned long idx, unsigned long fragnum)
{
  unsigned long *fragments;
  assert(chain);
  assert(idx < gt_array_size(chain->fragments));
  fragments = gt_array_get_space(chain->fragments);
  fragments[idx] = fragnum;
}

unsigned long chain_get_fragnum(const Chain *chain, unsigned long idx)
{
  assert(chain);
  assert(idx < gt_array_size(chain->fragments));
  return *(unsigned long*) gt_array_get(chain->fragments, idx);
}

unsigned long chain_size(const Chain *chain)
{
  assert(chain);
  return gt_array_size(chain->fragments);
}

void chain_delete(Chain *chain)
{
  if (!chain) return;
  gt_array_delete(chain->fragments);
  gt_free(chain);
}
