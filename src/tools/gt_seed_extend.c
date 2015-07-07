/*
  Copyright (c) 2015 Jörg Winkler <joerg.winkler@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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
#include "core/encseq_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "match/diagbandseed.h"
#include "tools/gt_seed_extend.h"

typedef struct {
  unsigned int kmerlen;
  unsigned int diagbandw;
  unsigned int mincoverage;
  unsigned int maxfreq;
  unsigned int correlation;
  unsigned int minalilen;
  unsigned int maxfrontdist;
  unsigned int minquality;
  bool mirror;
  bool verify;
  bool benchmark;
} GtSeedExtendArguments;

static void* gt_seed_extend_arguments_new(void)
{
  GtSeedExtendArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  return arguments;
}

static void gt_seed_extend_arguments_delete(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_free(arguments);
  }
}

static GtOptionParser* gt_seed_extend_option_parser_new(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] encseq_basename",
                            "Calculate local alignments using the seed and "
                            "extend algorithm.");

  /* -kmerlen */
  option = gt_option_new_uint_min("kmerlen", "k-mer length",
                                  &arguments->kmerlen, 14, 2);
  gt_option_parser_add_option(op, option);

  /* -diagbandw */
  option = gt_option_new_uint("diagbandw", "diagonal band width",
                              &arguments->diagbandw, 6);
  gt_option_parser_add_option(op, option);

  /* -mincoverage */
  option = gt_option_new_uint("mincoverage", "minimum coverage in two diagonal "
                              "bands", &arguments->mincoverage, 35);
  gt_option_parser_add_option(op, option);

  /* -maxfreq */
  option = gt_option_new_uint_min("maxfreq", "maximum frequency of a k-mer",
                                  &arguments->maxfreq, UINT_MAX, 1);
  gt_option_parser_add_option(op, option);

  /* -correlation */
  option = gt_option_new_uint("correlation", "percent similarity of the reads",
                              &arguments->correlation, 70);
  gt_option_parser_add_option(op, option);

  /* -minalilen */
  option = gt_option_new_uint_min("minalilen", "minimum alignment length",
                                  &arguments->minalilen, 1000, 1);
  gt_option_parser_add_option(op, option);

  /* -maxfrontdist */
  option = gt_option_new_uint("maxfrontdist", "maximum distance between 2 "
                              "fronts", &arguments->maxfrontdist, UINT_MAX);
  gt_option_parser_add_option(op, option);

  /* -minquality */
  option = gt_option_new_uint("minquality", "percent minimum alignment quality",
                              &arguments->minquality, 50);
  gt_option_parser_add_option(op, option);

  /* -mirror */
  option = gt_option_new_bool("mirror", "add reverse complement reads",
                              &arguments->mirror, false);
  gt_option_parser_add_option(op, option);

  /* -verify */
  option = gt_option_new_bool("verify", "check that k-mer seeds occur in the "
                              "sequences", &arguments->verify, false);
  gt_option_parser_add_option(op, option);

  /* -benchmark */
  option = gt_option_new_bool("benchmark", "measure time of different steps",
                              &arguments->benchmark, false);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_seed_extend_arguments_check(GT_UNUSED int rest_argc,
                                          void *tool_arguments, GtError *err)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (arguments->correlation > 100) {
    gt_error_set(err, "correlation must be <= 100");
    had_err = -1;
  }

  if (arguments->minquality > 100) {
    gt_error_set(err, "minquality must be <= 100");
    had_err = -1;
  }

  if (rest_argc == 1 && arguments->maxfreq == 1) {
    gt_error_set(err, "for 1 input file maxfreq must be >= 2 to find matching "
                 "k-mers");
    had_err = -1;
  }

  if (rest_argc > 2) {
    gt_error_set(err, "too many arguments (-help shows correct usage)");
    had_err = -1;
  } else if (rest_argc < 1) {
    gt_error_set(err, "at least one encseq index name must be specified "
      "(-help shows correct usage)");
    had_err = -1;
  }
  return had_err;
}

static int gt_seed_extend_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtEncseqLoader *encseq_loader;
  GtEncseq *aencseq, *bencseq;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);
  gt_assert(argc - parsed_args >= 1);

  /* load encseq A */
  encseq_loader = gt_encseq_loader_new();
  if (arguments->mirror)
    gt_encseq_loader_mirror(encseq_loader);
  gt_encseq_loader_enable_autosupport(encseq_loader);
  aencseq = gt_encseq_loader_load(encseq_loader, argv[parsed_args], err);
  if (!aencseq)
    had_err = -1;

  /* if there is a 2nd read set: load encseq B */
  if (!had_err) {
    if (argc - parsed_args == 2)
      bencseq = gt_encseq_loader_load(encseq_loader, argv[parsed_args+1], err);
    else
      bencseq = gt_encseq_ref(aencseq);
    if (!bencseq) {
      had_err = -1;
      gt_encseq_delete(aencseq);
    }
  }
  gt_encseq_loader_delete(encseq_loader);

  /* start algorithm */
  if (!had_err) {
    gt_seed_extend_run(aencseq, bencseq, (GtSeedExtend *)arguments);
    gt_encseq_delete(aencseq);
    gt_encseq_delete(bencseq);
  }
  return had_err;
}

GtTool* gt_seed_extend(void)
{
  return gt_tool_new(gt_seed_extend_arguments_new,
                     gt_seed_extend_arguments_delete,
                     gt_seed_extend_option_parser_new,
                     gt_seed_extend_arguments_check,
                     gt_seed_extend_runner);
}
