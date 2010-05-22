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

#include "core/ma.h"
#include "core/encseq_api.h"
#include "core/fasta_separator.h"
#include "core/str_array_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "tools/gt_seqdecode.h"

static GtOptionParser* gt_seqdecode_option_parser_new(GT_UNUSED
                                                      void *tool_arguments)
{
  GtOptionParser *op;

  /* init */
  op = gt_option_parser_new("sequence_file", "Decode encoded sequence_file.");

  gt_option_parser_set_min_max_args(op, 1, 1);

  return op;
}

static void output_sequence(GtEncseq *encseq)
{
  unsigned long i, j;
  gt_assert(encseq);
  for (i = 0; i < gt_encseq_num_of_sequences(encseq); i++) {
    GtEncseqReader *encseq_reader;
    unsigned long desclen;
    const char *desc;
    /* output description */
    gt_xfputc(GT_FASTA_SEPARATOR, stdout);
    desc = gt_encseq_description(encseq, &desclen, i);
    gt_xfwrite(desc, 1, desclen, stdout);
    gt_xfputc('\n', stdout);
    /* output sequence */
    encseq_reader = gt_encseq_create_reader_with_direction(encseq, true,
                                              gt_encseq_seqstartpos(encseq, i));
    /* XXX: make this more efficient by writing in a buffer first and the
       showing the result */
    for (j = 0; j < gt_encseq_seqlength(encseq, i);  j++) {
      gt_xfputc(gt_encseq_reader_next_decoded_char(encseq_reader), stdout);
    }
    gt_xfputc('\n', stdout);
    gt_encseq_reader_delete(encseq_reader);
  }
}

static int decode_sequence_file(const char *seqfile, GtError *err)
{
  GtEncseqLoader *encseq_loader;
  GtEncseq *encseq;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(seqfile);
  encseq_loader = gt_encseq_loader_new();
  if (!(encseq = gt_encseq_loader_load(encseq_loader, seqfile, err)))
    had_err = -1;
  if (!had_err)
    output_sequence(encseq);
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(encseq_loader);
  return had_err;
}

static int gt_seqdecode_runner(GT_UNUSED int argc, const char **argv,
                               int parsed_args, GT_UNUSED void *tool_arguments,
                               GtError *err)
{
  gt_error_check(err);
  return decode_sequence_file(argv[parsed_args], err);
}

GtTool* gt_seqdecode(void)
{
  return gt_tool_new(NULL,
                     NULL,
                     gt_seqdecode_option_parser_new,
                     NULL,
                     gt_seqdecode_runner);
}
