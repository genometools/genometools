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

struct GtChain {
  GtArray *fragments;
  long score;
};

GtChain* gt_chain_new(void)
{
  GtChain *chain;
  chain = gt_calloc((size_t) 1, sizeof *chain);
  chain->fragments = gt_array_new(sizeof (unsigned long));
  return chain;
}

void gt_chain_reset(GtChain *chain)
{
  gt_assert(chain);
  chain->score = 0;
  gt_array_reset(chain->fragments);
}

long gt_chain_get_score(const GtChain *chain)
{
  gt_assert(chain);
  return chain->score;
}

void gt_chain_set_score(GtChain *chain, long score)
{
  gt_assert(chain);
  chain->score = score;
}

void gt_chain_add_fragnum(GtChain *chain, unsigned long fragnum)
{
  gt_assert(chain);
  gt_array_add(chain->fragments, fragnum);
}

void gt_chain_set_fragnum(GtChain *chain, unsigned long idx,
                          unsigned long fragnum)
{
  unsigned long *fragments;
  gt_assert(chain);
  gt_assert(idx < gt_array_size(chain->fragments));
  fragments = gt_array_get_space(chain->fragments);
  fragments[idx] = fragnum;
}

unsigned long gt_chain_get_fragnum(const GtChain *chain, unsigned long idx)
{
  gt_assert(chain);
  gt_assert(idx < gt_array_size(chain->fragments));
  return *(unsigned long*) gt_array_get(chain->fragments, idx);
}

unsigned long gt_chain_size(const GtChain *chain)
{
  gt_assert(chain);
  return gt_array_size(chain->fragments);
}

void gt_chain_delete(GtChain *chain)
{
  if (!chain) return;
  gt_array_delete(chain->fragments);
  gt_free(chain);
}
