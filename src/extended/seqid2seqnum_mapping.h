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

#ifndef SEQID2SEQNUM_MAPPING_H
#define SEQID2SEQNUM_MAPPING_H

#include "core/bioseq.h"
#include "core/range_api.h"
#include "core/seq_col.h"

/* Helper class which maps a sequence id to a sequence number. */
typedef struct GtSeqid2SeqnumMapping GtSeqid2SeqnumMapping;

/* Create a new <GtSeqid2SeqnumMapping*> object which is able to map the
   sequence ids given in <seqcol> to the corresponding sequence and file
   numbers.
   The descriptions are parsed for ``description ranges'' and the corresponding
   offsets are stored (descriptions without a description range have an offset
   of 1). */
GtSeqid2SeqnumMapping* gt_seqid2seqnum_mapping_new_seqcol(GtSeqCol *seq_col,
                                                          GtError*);
/* Create a new <GtSeqid2SeqnumMapping*> object which is able to map the
   sequence ids given in <bioseq> to the corresponding sequence numbers.
   The descriptions are parsed for ``description ranges'' and the corresponding
   offsets are stored (descriptions without a description range have an offset
   of 1). */
GtSeqid2SeqnumMapping* gt_seqid2seqnum_mapping_new_bioseq(GtBioseq *bioseq,
                                                          GtError*);
/* Delete the given <seqid2seqnum_mapping> object. */
void                   gt_seqid2seqnum_mapping_delete(GtSeqid2SeqnumMapping
                                                      *seqid2seqnum_mapping);
/* Use <seqid2seqnum_mapping> to map the given <seqid> (in the given <range>)
   and store the result in <seqnum>. The corresponding offset is stored in
   <offset>. If the mapping is not possible, -1 is returned and <err> is set
   accordingly. */
int                    gt_seqid2seqnum_mapping_map(GtSeqid2SeqnumMapping
                                                   *seqid2seqnum_mapping,
                                                   const char *seqid,
                                                   const GtRange *range,
                                                   unsigned long *seqnum,
                                                   unsigned long *filenum,
                                                   unsigned long *offset,
                                                   GtError *err);

#endif
