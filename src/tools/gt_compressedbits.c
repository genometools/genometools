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
#include <stdio.h>

#include "core/fa.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/str_api.h"
#include "core/str_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/compressed_bitsequence.h"
#include "extended/popcount_tab.h"
#include "tools/gt_compressedbits.h"

typedef struct {
  unsigned int samplerate;
  GtUword size,
                benches;
  bool fill_random,
       check_consistency;
  GtStr *filename;
  GtOption *size_op,
           *filename_op,
           *rand_op;
} GtCompressdbitsArguments;

static void* gt_compressedbits_arguments_new(void)
{
  GtCompressdbitsArguments *arguments =
    gt_calloc((size_t) 1, sizeof *arguments);
  arguments->filename = gt_str_new();
  return arguments;
}

static void gt_compressedbits_arguments_delete(void *tool_arguments)
{
  GtCompressdbitsArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->filename);
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
  op = gt_option_parser_new("[option ...]",
                       "Testing compressed bitsequence, save to disk, reload.");

  /* -size */
  option = gt_option_new_uword("size",
                               "size of GtBitsequence to create "
                               "(words 32/64 bit)",
                               &arguments->size, 20UL);
  gt_option_parser_add_option(op, option);
  arguments->size_op = gt_option_ref(option);

  /* -samplerate */
  option = gt_option_new_uint("samplerate",
                              "samplerate of random GtBitsequence to test",
                              &arguments->samplerate, 32U);
  gt_option_parser_add_option(op, option);

  /* -rand */
  option = gt_option_new_bool("rand", "create random bitvector",
                              &arguments->fill_random, false);
  gt_option_parser_add_option(op, option);
  arguments->rand_op = gt_option_ref(option);

  /* -check */
  option = gt_option_new_bool("check", "compare original with compressed and "
                              "loaded from file",
                              &arguments->check_consistency, false);
  gt_option_parser_add_option(op, option);
  arguments->rand_op = gt_option_ref(option);

  /* -input */
  option = gt_option_new_filename(
                                "input",
                                "load vector from file, format is as follows:\n"
                                "[ULL size in bits][[ULL bits]...]\n"
                                " not usable with -size and -rand",
                                arguments->filename);
  gt_option_parser_add_option(op, option);
  arguments->filename_op = gt_option_ref(option);
  gt_option_exclude(arguments->filename_op, arguments->size_op);
  gt_option_exclude(arguments->filename_op, arguments->rand_op);
  /* -benches */
  option = gt_option_new_uword("benches",
                               "number of function calls to benchmark",
                               &arguments->benches, 100000UL);
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
  GtUword idx;
  GtUint64 num_of_bits = 0ULL;
  GtBitsequence *bits = NULL;
  GtCompressedBitsequence *cbs = NULL, *read_cbs = NULL;
  GtStr *filename = gt_str_new();
  FILE *fp = NULL;

  gt_error_check(err);
  gt_assert(arguments);
  gt_assert(argc == parsed_args);

  if (gt_option_is_set(arguments->filename_op)) {
    FILE *file = NULL;
    gt_assert(arguments->filename != NULL);

    file = gt_xfopen(gt_str_get(arguments->filename), "r");
    if ((size_t) 1 != gt_xfread(&num_of_bits,
                                sizeof (num_of_bits), (size_t) 1, file)) {
      had_err = -1;
    }
    if (!had_err) {
      gt_log_log("bits to read: "GT_LLU"", num_of_bits);
      arguments->size = (GtUword) GT_NUMOFINTSFORBITS(num_of_bits);
      bits = gt_malloc(sizeof (*bits) * arguments->size);
      if ((size_t) arguments->size !=
          gt_xfread(bits, sizeof (*bits),
                    (size_t) arguments->size, file)) {
        had_err = -1;
      }
    }
    gt_xfclose(file);
  }
  else {
    bits = gt_calloc(sizeof (*bits), (size_t) arguments->size);
    num_of_bits = (GtUint64) (GT_INTWORDSIZE * arguments->size);

    if (arguments->fill_random) {
      for (idx = 0; idx < arguments->size; idx++) {
        bits[idx] =
          (GtBitsequence) (0xAAAAAAAAAAAAAAAAULL ^ gt_rand_max(ULONG_MAX));
      }
    }
    else {
      for (idx = 0; idx < arguments->size; idx++)
        bits[idx] = (GtBitsequence) (0xAAAAAAAAAAAAAAAAULL ^ idx);
    }
  }

  if (!had_err) {
    fp = gt_xtmpfp(filename);
    gt_fa_xfclose(fp);
    fp = NULL;

    gt_log_log("filename: %s", gt_str_get(filename));
    gt_log_log("size in words: "GT_WU"", arguments->size);
    cbs = gt_compressed_bitsequence_new(
                            bits, arguments->samplerate,
                            (GtUword) num_of_bits);
    gt_log_log("original size in MB: %2.3f",
               (sizeof (*bits) * arguments->size) / (1024.0 * 1024.0));
    gt_log_log("compressed size in MB: %2.3f",
               gt_compressed_bitsequence_size(cbs) / (1024.0 * 1024.0));
    gt_log_log("popcount table size thereof in MB: %2.3f",
               gt_popcount_tab_calculate_size(15U) / (1024.0 * 1024.0));
    had_err = gt_compressed_bitsequence_write(cbs, gt_str_get(filename), err);
  }
  if (!had_err)
  {
    read_cbs =
      gt_compressed_bitsequence_new_from_file(gt_str_get(filename), err);
    if (read_cbs == NULL)
      had_err = -1;
  }
  if (!had_err && bits != NULL && arguments->check_consistency) {
    for (idx = 0; (GtUint64) idx < num_of_bits; ++idx) {
      int GT_UNUSED bit = gt_compressed_bitsequence_access(read_cbs, idx);
      int GT_UNUSED original = GT_ISIBITSET(bits, idx) ? 1 : 0;
      gt_assert(gt_compressed_bitsequence_access(cbs, idx) == bit);
      gt_assert(original == bit);
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
