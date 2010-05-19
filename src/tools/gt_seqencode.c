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
#include "core/str_array_api.h"
#include "core/unused_api.h"
#include "tools/gt_seqencode.h"

static GtOptionParser* gt_seqencode_option_parser_new(GT_UNUSED
                                                      void *tool_arguments)
{
  GtOptionParser *op;

  /* init */
  op = gt_option_parser_new("sequence_file",
                            "Encode sequence_file efficiently.");

  gt_option_parser_set_min_max_args(op, 1, 1);

  return op;
}

static int encode_sequence_file(const char *seqfile, GtError *err)
{
  GtEncseqEncoder *encseq_encoder;
  GtStrArray *seqfiles;
  int had_err;
  gt_error_check(err);
  gt_assert(seqfile);
  encseq_encoder = gt_encseq_encoder_new();
  seqfiles = gt_str_array_new();
  gt_str_array_add_cstr(seqfiles, seqfile);
  had_err = gt_encseq_encoder_encode(encseq_encoder, seqfiles, seqfile, err);
  gt_str_array_delete(seqfiles);
  gt_encseq_encoder_delete(encseq_encoder);
  return had_err;
}

static int gt_seqencode_runner(GT_UNUSED int argc, const char **argv,
                               int parsed_args, GT_UNUSED void *tool_arguments,
                               GtError *err)
{
  gt_error_check(err);
  return encode_sequence_file(argv[parsed_args], err);
}

GtTool* gt_seqencode(void)
{
  return gt_tool_new(NULL,
                     NULL,
                     gt_seqencode_option_parser_new,
                     NULL,
                     gt_seqencode_runner);
}
