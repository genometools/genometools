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
#include "core/unused_api.h"

#include "match/shu-encseq-gc.h"

static inline void calculate_gc(const GtEncseq *encseq,
                                double *gc_content,
                                bool with_special,
                                unsigned long seq_idx,
                                unsigned long gc_count,
                                unsigned long at_count)
{
  if (with_special)
  {
    gt_assert(seq_idx < gt_encseq_num_of_sequences(encseq));
    gc_content[seq_idx] =
      (double) gc_count / (double) gt_encseq_seqlength(encseq, seq_idx);
  }
  else
  {
    gc_content[seq_idx] = (double) gc_count / (double) (gc_count + at_count);
  }
}

double *gt_encseq_get_gc(const GtEncseq *encseq,
                         bool with_special,
                         bool calculate,
                         GT_UNUSED GtError *err)
{
  GtEncseqReader *reader;
  GtAlphabet *alphabet;
  double *gc_content;
  /* unit = file or sequence depending on per_file */
  unsigned long char_idx, totallength, max_unit,
                seq_idx = 0,
                nextsep = 0,
                at_count = 0,
                gc_count = 0,
                default_count = 0;
  bool is_mirrored_encseq;
  GtUchar acgt[8], current_c;

  alphabet = gt_encseq_alphabet(encseq);
  gt_assert(gt_alphabet_is_dna(alphabet));
  gt_alphabet_encode_seq(alphabet, acgt,
                         "aAtTcCgG", 8UL);
  totallength = gt_encseq_total_length(encseq);
  reader = gt_encseq_create_reader_with_readmode(encseq,
                                                 GT_READMODE_FORWARD,
                                                 0);
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

  nextsep = gt_encseq_seqstartpos(encseq, seq_idx) +
            gt_encseq_seqlength(encseq, seq_idx);

  for (char_idx = 0; char_idx < totallength; char_idx++)
  {
    if (nextsep == char_idx)
    {
      if (calculate)
      {
        calculate_gc(encseq,
                     gc_content,
                     with_special,
                     seq_idx,
                     gc_count,
                     at_count);
      }
      else
      {
        gc_content[seq_idx] = (double) gc_count;
      }

      seq_idx++;

      nextsep = gt_encseq_seqstartpos(encseq, seq_idx) +
                gt_encseq_seqlength(encseq, seq_idx);

      gt_encseq_reader_reinit_with_readmode(reader,
                                            encseq,
                                            GT_READMODE_FORWARD,
                                            char_idx + 1UL);
      gc_count = at_count = default_count = 0UL;
      continue;
    }
    current_c = gt_encseq_reader_next_encoded_char(reader);
    if (current_c == acgt[0] ||
        current_c == acgt[1] ||
        current_c == acgt[2] ||
        current_c == acgt[3])
    {
       at_count++;
    }
    else
    {
      if (current_c == acgt[4] ||
          current_c == acgt[5] ||
          current_c == acgt[6] ||
          current_c == acgt[7])
      {
         gc_count++;
      }
      else
      {
        default_count++;
      }
    }
  }
  if (calculate)
  {
    calculate_gc(encseq,
                 gc_content,
                 with_special,
                 seq_idx,
                 gc_count,
                 at_count);
  }
  else
  {
    gc_content[seq_idx] = (double) gc_count;
  }
  gt_encseq_reader_delete(reader);
  if (is_mirrored_encseq)
  {
    unsigned long double_max_unit = GT_MULT2(max_unit);
    for (seq_idx = 0; seq_idx < max_unit; seq_idx++)
    {
      gc_content[double_max_unit - seq_idx - 1] =
        gc_content[seq_idx];
    }
  }
  return gc_content;
}

int gt_encseq_gc_unit_test(GtError *err)
{
  int had_err = 0;
  double *results;
  GtEncseqBuilder *eb;
  GtEncseq *encseq;
  const char testseq1[] = "aaaaaa",
             testseq2[] = "cccccc",
             testseq3[] = "acgtacgt",
             testseq4[] = "acgtn";
          /* testseq5[] = "xxxxn"; */
  GtAlphabet *alpha;
/*GtError *tmp_err;*/

  gt_error_check(err);

  alpha = gt_alphabet_new_dna();

  /* test a-seq */
  eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_create_ssp_tab(eb);
  gt_encseq_builder_enable_description_support(eb);
  gt_encseq_builder_add_cstr(eb, testseq1, 6UL, "only a");
  encseq = gt_encseq_builder_build(eb, err);
  if ((results = gt_encseq_get_gc(encseq,
                                  false,
                                  true,
                                  err)) != NULL)
  {
    gt_ensure(had_err, gt_double_equals_double(results[0], 0.0));
  }
  else
  {
    had_err = -1;
  }
  gt_free(results);
  gt_encseq_builder_delete(eb);
  gt_encseq_delete(encseq);

  if (!had_err)
  {
    /* test c-seq */
    eb = gt_encseq_builder_new(alpha);
    gt_encseq_builder_create_ssp_tab(eb);
    gt_encseq_builder_enable_description_support(eb);
    gt_encseq_builder_add_cstr(eb, testseq2, 6UL, "only c");
    encseq = gt_encseq_builder_build(eb, err);
    if ((results = gt_encseq_get_gc(encseq,
                                    false,
                                    true,
                                    err)) != NULL)
    {
      gt_ensure(had_err, gt_double_equals_one(results[0]));
    }
    else
    {
      had_err = -1;
    }
    gt_free(results);
    gt_encseq_builder_delete(eb);
    gt_encseq_delete(encseq);
  }

  if (!had_err)
  {
    /* test dna-seq and dna+special-seq*/
    eb = gt_encseq_builder_new(alpha);
    gt_encseq_builder_create_ssp_tab(eb);
    gt_encseq_builder_enable_description_support(eb);
    gt_encseq_builder_add_cstr(eb, testseq3, 8UL, "0.5");
    gt_encseq_builder_add_cstr(eb, testseq4, 5UL, "0.5+special");
    encseq = gt_encseq_builder_build(eb, err);
    if ((results = gt_encseq_get_gc(encseq,
                                    false,
                                    true,
                                    err)) != NULL)
    {
      gt_ensure(had_err, gt_double_equals_double(results[0], 0.5));
      gt_ensure(had_err, gt_double_equals_double(results[1], 0.5));
    }
    else
    {
      had_err = -1;
    }
    gt_free(results);

    if (!had_err)
    {
      /* count special chars */
      if ((results = gt_encseq_get_gc(encseq,
                                      true,
                                      true,
                                      err)) != NULL)
      {
        gt_ensure(had_err, gt_double_equals_double(results[0], 0.5));
        gt_ensure(had_err, gt_double_equals_double(results[1], (2.0/5.0)));
      }
      else
      {
        had_err = -1;
      }
      gt_free(results);
    }

    gt_encseq_builder_delete(eb);
    gt_encseq_delete(encseq);
  }

  if (!had_err)
    {
      /* test dna-seq and dna+special-seq*/
      eb = gt_encseq_builder_new(alpha);
      gt_encseq_builder_create_ssp_tab(eb);
      gt_encseq_builder_enable_description_support(eb);
      gt_encseq_builder_add_cstr(eb, testseq3, 8UL, "0.5");
      gt_encseq_builder_add_cstr(eb, testseq4, 5UL, "0.5+special");
      encseq = gt_encseq_builder_build(eb, err);
      /*add mirrored sequence*/
      had_err = gt_encseq_mirror(encseq, err);
      /* sequence wise */
      if ((results = gt_encseq_get_gc(encseq,
                                      false,
                                      true,
                                      err)) != NULL)
      {
        gt_ensure(had_err, gt_double_equals_double(results[0], 0.5));
        gt_ensure(had_err, gt_double_equals_double(results[1], 0.5));
        gt_ensure(had_err, gt_double_equals_double(results[2], 0.5));
        gt_ensure(had_err, gt_double_equals_double(results[3], 0.5));
      }
      else
      {
        had_err = -1;
      }
      gt_free(results);

      if (!had_err)
      {
        /* count special chars */
        if ((results = gt_encseq_get_gc(encseq,
                                        true,
                                        true,
                                        err)) != NULL)
        {
          gt_ensure(had_err, gt_double_equals_double(results[0], 0.5));
          gt_ensure(had_err, gt_double_equals_double(results[1], (2.0/5.0)));
          gt_ensure(had_err, gt_double_equals_double(results[2], (2.0/5.0)));
          gt_ensure(had_err, gt_double_equals_double(results[3], 0.5));
        }
        else
        {
          had_err = -1;
        }
        gt_free(results);
      }
      gt_encseq_builder_delete(eb);
      gt_encseq_delete(encseq);
    }
    gt_alphabet_delete(alpha);
  return had_err;
}
