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

#include "core/ma_api.h"
#include "gth/chain_collection.h"

struct GthChainCollection {
  GtArray *chains;
};

GthChainCollection* gth_chain_collection_new(void)
{
  GthChainCollection *chain_collection = gt_malloc(sizeof *chain_collection);
  chain_collection->chains = gt_array_new(sizeof (GthChain*));
  return chain_collection;
}

void gth_chain_collection_delete(GthChainCollection *chain_collection)
{
  unsigned long i;
  if (!chain_collection) return;
  for (i = 0; i < gt_array_size(chain_collection->chains); i++)
    gth_chain_delete(*(GthChain**) gt_array_get(chain_collection->chains, i));
  gt_array_delete(chain_collection->chains);
  gt_free(chain_collection);
}

void gth_chain_collection_add(GthChainCollection *chain_collection,
                              GthChain *chain)
{
  gt_assert(chain_collection && chain);
  gt_array_add(chain_collection->chains, chain);
}

static int compare_chains(const void *dataA, const void *dataB)
{
  GthChain *chainA = *(GthChain**) dataA;
  GthChain *chainB = *(GthChain**) dataB;

  gt_assert(chainA->ref_file_num == chainB->ref_file_num);

  /* sort after genomic file number */
  if (chainA->gen_file_num < chainB->gen_file_num)
    return -1;
  if (chainA->gen_file_num > chainB->gen_file_num)
    return 1;

  /* genomic file number is the same,
     sort after genomic sequence number */
  if (chainA->gen_seq_num < chainB->gen_seq_num)
    return -1;
  if (chainA->gen_seq_num > chainB->gen_seq_num)
    return 1;

  /* genomic sequence number is also the same,
     sort after reference sequence number */
  if (chainA->ref_seq_num < chainB->ref_seq_num)
    return -1;
  if (chainA->ref_seq_num > chainB->ref_seq_num)
    return 1;

  /* reference sequence number is also the same,
     sort after reference sequence coverage */
  if (chainA->refseqcoverage < chainB->refseqcoverage)
    return -1;
  if (chainA->refseqcoverage > chainB->refseqcoverage)
    return 1;

  return 0;
}

void gth_chain_collection_sort(GthChainCollection *chain_collection)
{
  gt_assert(chain_collection);
  if (gt_array_size(chain_collection->chains)) {
    qsort(gt_array_get_space(chain_collection->chains),
          gt_array_size(chain_collection->chains), sizeof (GthChain*),
          compare_chains);
  }
}

unsigned long gth_chain_collection_size(const GthChainCollection
                                        *chain_collection)
{
  gt_assert(chain_collection);
  return gt_array_size(chain_collection->chains);
}

GthChain* gth_chain_collection_get(const GthChainCollection *chain_collection,
                                   unsigned long idx)
{
  gt_assert(chain_collection);
  return *(GthChain**) gt_array_get(chain_collection->chains, idx);
}
