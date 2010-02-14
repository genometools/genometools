/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/hashmap_api.h"
#include "core/ma_api.h"
#include "extended/seqid2seqnum_mapping.h"

struct GtSeqid2SeqnumMapping {
  GtHashmap *map;
  unsigned long *seqnums,
                *offsets;
};

static int fill_mapping(GtSeqid2SeqnumMapping *mapping, GtBioseq *bioseq,
                        GtError *err)
{
  gt_error_check(err);
  gt_assert(mapping && bioseq);

  /* XXX */

  return 0;
}

GtSeqid2SeqnumMapping* gt_seqid2seqnum_mapping_new(GtBioseq *bioseq,
                                                   GtError *err)
{
  unsigned long num_of_seqs;
  GtSeqid2SeqnumMapping *mapping;
  gt_error_check(err);
  gt_assert(bioseq);
  mapping = gt_malloc(sizeof *mapping);
  mapping->map = gt_hashmap_new(GT_HASH_DIRECT, gt_free_func, NULL);
  num_of_seqs = gt_bioseq_number_of_sequences(bioseq);
  mapping->seqnums = gt_malloc(num_of_seqs * sizeof *mapping->seqnums);
  mapping->offsets = gt_malloc(num_of_seqs * sizeof *mapping->offsets);
  if (fill_mapping(mapping, bioseq, err)) {
    gt_seqid2seqnum_mapping_delete(mapping);
    return NULL;
  }
  return mapping;
}

void gt_seqid2seqnum_mapping_delete(GtSeqid2SeqnumMapping *mapping)
{
  if (!mapping) return;
  gt_free(mapping->offsets);
  gt_free(mapping->seqnums);
  gt_hashmap_delete(mapping->map);
  gt_free(mapping);
}

int gt_seqid2seqnum_mapping_map(GtSeqid2SeqnumMapping *mapping,
                                const char *seqid, unsigned long *seqnum,
                                unsigned long *offset, GtError *err)
{
  gt_error_check(err);
  gt_assert(mapping && seqid && seqnum && offset);

  /* XXX */
  *offset = 1;

  return 0;
}
