/*
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include <stdbool.h>

#include "core/alphabet_api.h"
#include "core/divmodmul.h"
#include "core/encseq_api.h"
#include "core/ensure.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/safearith.h"
#include "core/unused_api.h"

#include "match/shu-encseq-gc.h"

static inline void calculate_rel_gc(const GtEncseq *encseq,
                                    double *gc_content,
                                    unsigned long seq_idx,
                                    unsigned long gc_count)
{
  unsigned long length;
  gt_assert(seq_idx < gt_encseq_num_of_sequences(encseq));
  length = gt_encseq_seqlength(encseq, seq_idx);
  gt_assert(gc_count <= length);
  gc_content[seq_idx] = (double) gc_count / (double) length;
  gt_assert(gt_double_compare(gc_content[seq_idx], 0.0) != -1);
  gt_assert(gt_double_compare(gc_content[seq_idx], 1.0) != 1);
}

unsigned long *gt_encseq_gc_count(const GtEncseq *encseq)
{
  GtEncseqReader *reader;
  GtAlphabet *alphabet;
  unsigned long *gc_count_per_seq,
                char_idx, totallength, max_unit,
                seq_idx = 0,
                nextsep = 0;
  bool is_mirrored_encseq;
  GtUchar cCgG_encoded[4], current_c;

  alphabet = gt_encseq_alphabet(encseq);
  gt_assert(gt_alphabet_is_dna(alphabet));
  gt_alphabet_encode_seq(alphabet, cCgG_encoded, "cCgG", 4UL);
  totallength = gt_encseq_total_length(encseq);
  reader = gt_encseq_create_reader_with_readmode(encseq,
                                                 GT_READMODE_FORWARD,
                                                 0);
  is_mirrored_encseq = gt_encseq_is_mirrored(encseq);
  if (is_mirrored_encseq) {
    max_unit = GT_DIV2(gt_encseq_num_of_sequences(encseq));
    gc_count_per_seq = gt_calloc((size_t) GT_MULT2(max_unit),
                                 sizeof (*gc_count_per_seq));
  }
  else {
    max_unit = gt_encseq_num_of_sequences(encseq);
    gc_count_per_seq = gt_calloc((size_t) max_unit, sizeof (*gc_count_per_seq));
  }

  nextsep = gt_encseq_seqstartpos(encseq, seq_idx) +
            gt_encseq_seqlength(encseq, seq_idx);

  for (char_idx = 0; char_idx < totallength; char_idx++)
  {
    if (nextsep == char_idx)
    {
      gt_assert(gc_count_per_seq[seq_idx] <=
                gt_encseq_seqlength(encseq, seq_idx));
      seq_idx++;

      nextsep = gt_encseq_seqstartpos(encseq, seq_idx) +
                gt_encseq_seqlength(encseq, seq_idx);

      gt_encseq_reader_reinit_with_readmode(reader,
                                            encseq,
                                            GT_READMODE_FORWARD,
                                            char_idx + 1UL);
      char_idx++;
    }
    current_c = gt_encseq_reader_next_encoded_char(reader);
    if (current_c == cCgG_encoded[0] || current_c == cCgG_encoded[1] ||
        current_c == cCgG_encoded[2] || current_c == cCgG_encoded[3]) {
       gc_count_per_seq[seq_idx]++;
       gt_assert(gc_count_per_seq[seq_idx] != 0);
    }
  }
  gt_encseq_reader_delete(reader);
  return gc_count_per_seq;
}

double *gt_encseq_get_rel_gc(const GtEncseq *encseq,
                             GT_UNUSED GtError *err)
{
  double *gc_content;
  unsigned long *gc_count_per_seq,
                max_unit,
                seq_idx = 0;
  bool is_mirrored_encseq;

  is_mirrored_encseq = gt_encseq_is_mirrored(encseq);
  if (is_mirrored_encseq)
  {
    max_unit = GT_DIV2(gt_encseq_num_of_sequences(encseq));
    gc_content = gt_calloc((size_t) GT_MULT2(max_unit), sizeof (double));
  }
  else
  {
    max_unit = gt_encseq_num_of_sequences(encseq);
    gc_content = gt_calloc((size_t) max_unit, sizeof (double));
  }

  gc_count_per_seq = gt_encseq_gc_count(encseq);
  for (seq_idx = 0; seq_idx < max_unit; seq_idx++) {
    calculate_rel_gc(encseq,
                     gc_content,
                     seq_idx,
                     gc_count_per_seq[seq_idx]);
  }
  if (is_mirrored_encseq)
  {
    unsigned long double_max_unit = GT_MULT2(max_unit);
    for (seq_idx = 0; seq_idx < max_unit; seq_idx++)
    {
      gc_content[double_max_unit - seq_idx - 1] = gc_content[seq_idx];
    }
  }
  return gc_content;
}

int gt_encseq_gc_unit_test(GT_UNUSED GtError *err)
{
  int had_err = 0;
  /* XXX write new tests for the new functions */
  return had_err;
}
