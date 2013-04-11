/*
  Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>

#include "core/fa.h"
#include "core/intbits.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/str_api.h"
#include "core/unused_api.h"
#include "extended/compressed_bitsequence.h"
#include "tools/gt_compressedbits.h"

typedef struct {
  unsigned int size;
} GtCompressdbitsArguments;

static void* gt_compressedbits_arguments_new(void)
{
  GtCompressdbitsArguments *arguments =
    gt_calloc((size_t) 1, sizeof *arguments);
  return arguments;
}

static void gt_compressedbits_arguments_delete(void *tool_arguments)
{
  GtCompressdbitsArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_free(arguments);
  }
}

static GtOptionParser* gt_compressedbits_option_parser_new(void *tool_arguments)
{
  GtCompressdbitsArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]", /* XXX */
                            "Testing compressed bitsequence."); /* XXX */

  /* -size */
  option = gt_option_new_uint("size", "size of random GtBitsequence to test",
                              &arguments->size, 20U);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_compressedbits_runner(GT_UNUSED int argc,
                                   GT_UNUSED const char **argv,
                                   GT_UNUSED int parsed_args,
                                   void *tool_arguments,
                                   GtError *err)
{
  GtCompressdbitsArguments *arguments = tool_arguments;
  int had_err = 0;
  unsigned int samplerate;
  unsigned long idx;
  GtBitsequence *bits = NULL;
  GtCompressedBitsequence *cbs, *read_cbs = NULL;
  GtStr *filename = gt_str_new();
  FILE *fp;

  gt_error_check(err);
  gt_assert(arguments);
  gt_assert(argc == parsed_args);

  bits = gt_malloc(sizeof (*bits) * arguments->size);

  for (idx = 0; idx < (unsigned long) arguments->size; idx++) {
    bits[idx] = (GtBitsequence) gt_rand_max(ULONG_MAX);
  }

  fp = gt_xtmpfp(filename);
  gt_fa_xfclose(fp);
  fp = NULL;

  fprintf(stderr, "filename: %s\n", gt_str_get(filename));
  samplerate = arguments->size / 100;
  if (samplerate == 0) {
    samplerate = 5U;
  }
  cbs =
    gt_compressed_bitsequence_new(
                            bits, samplerate,
                            (unsigned long) (GT_LOGWORDSIZE * arguments->size));
  had_err = gt_compressed_bitsequence_write(cbs, gt_str_get(filename), err);
  if (!had_err)
  {
    read_cbs =
      gt_compressed_bitsequence_new_from_file(gt_str_get(filename), err);
    if (read_cbs == NULL)
      had_err = -1;
  }
  if (!had_err) {
    unsigned long num_of_bits =
      (unsigned long) (GT_LOGWORDSIZE * arguments->size);
    for (idx = 0; idx < num_of_bits; ++idx) {
      gt_assert(gt_compressed_bitsequence_access(cbs, idx) ==
                gt_compressed_bitsequence_access(read_cbs, idx));
    }
  }
  gt_compressed_bitsequence_delete(cbs);
  gt_compressed_bitsequence_delete(read_cbs);
  gt_free(bits);
  gt_str_delete(filename);
  return had_err;
}

GtTool* gt_compressedbits(void)
{
  return gt_tool_new(gt_compressedbits_arguments_new,
                     gt_compressedbits_arguments_delete,
                     gt_compressedbits_option_parser_new,
                     NULL,
                     gt_compressedbits_runner);
}
