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

#include "core/basename_api.h"
#include "core/fa.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/logger_api.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/n_r_encseq.h"
#include "match/xdrop.h"
#include "tools/gt_condenser_compress.h"

typedef struct {
  GtEncseq              *input_es;
  GtStr                 *indexname;
  GtXdropArbitraryscores scores;
  GtUword                minalignlength,
                         max_kmer_poss,
                         initsize;
  GtWord                 xdrop;
  unsigned int           kmersize,
                         windowsize;
  bool                   opt,
                         verbose;
} GtCondenserCompressArguments;

static void* gt_condenser_compress_arguments_new(void)
{
  GtCondenserCompressArguments *arguments = gt_calloc((size_t) 1,
                                                      sizeof *arguments);
  arguments->indexname = gt_str_new();
  arguments->input_es = NULL;
  return arguments;
}

static void gt_condenser_compress_arguments_delete(void *tool_arguments)
{
  GtCondenserCompressArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->indexname);
    gt_encseq_delete(arguments->input_es);
    gt_free(arguments);
  }
}

#define GT_CONDENSER_STR(X) #X
#define GT_CONDENSER_XSTR(X) GT_CONDENSER_STR(X)

static GtOptionParser*
gt_condenser_compress_option_parser_new(void *tool_arguments)
{
  GtCondenserCompressArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[options] INPUTENCSEQ",
                            "Compresses a GtEncseq to a UniqueEncseq.");

  /* -indexname */
  option = gt_option_new_string("indexname",
                                "path and basename of files to store",
                                arguments->indexname, NULL);
  gt_option_parser_add_option(op, option);

  /* -kmersize */
  option = gt_option_new_uint_min("kmersize",
                                  "kmer-size used for the seeds, default "
                                  "depends on alphabet size",
                                  &arguments->kmersize, GT_UNDEF_UINT, 2U);
  gt_option_parser_add_option(op, option);

  /* -windowsize */
  option = gt_option_new_uint("windowsize",
                              "Size of window in which to search for hit pairs "
                              "of kmers, has to be larger than kmersize" ,
                              &arguments->windowsize, GT_UNDEF_UINT);
  gt_option_parser_add_option(op, option);

  /* -initsize */
  option = gt_option_new_uword("initsize",
                               "length of inital unique database in bases, "
                               "should be larger than -alignlength",
                               &arguments->initsize, GT_UNDEF_UWORD);
  gt_option_parser_add_option(op, option);

  /* -alignlength */
  option = gt_option_new_uword("alignlength",
                               "required minimal length of an xdrop-alignment, "
                               "should be larger than -windowsize",
                               &arguments->minalignlength, GT_UNDEF_UWORD);
  gt_option_parser_add_option(op, option);

  /* -maxkmerpos */
  option = gt_option_new_uword_min("maxkmerpos",
                                   "maximal number of positions stored per "
                                   "kmer, min "
                                   GT_CONDENSER_XSTR(GT_NRENCSEQ_MIN_KMER_POS),
                                   &arguments->max_kmer_poss,
                                   (GtUword) 100,
                                   (GtUword) GT_NRENCSEQ_MIN_KMER_POS);
  gt_option_parser_add_option(op, option);

  /* -mat */
  option = gt_option_new_int("mat",
                             "matchscore for extension-alignment, "
                             "requirements: mat > mis, mat > 2ins, mat > 2del",
                             &arguments->scores.mat, 2);
  gt_option_parser_add_option(op, option);

  /* -mis */
  option = gt_option_new_int("mis",
                             "mismatchscore for extension-alignment, ",
                             &arguments->scores.mis, -1);
  gt_option_parser_add_option(op, option);

  /* -ins */
  option = gt_option_new_int("ins",
                             "insertionscore for extension-alignment",
                             &arguments->scores.ins, -2);
  gt_option_parser_add_option(op, option);

  /* -del */
  option = gt_option_new_int("del",
                             "deletionscore for extension-alignment",
                             &arguments->scores.del, -2);
  gt_option_parser_add_option(op, option);

  /* -xdrop */
  option = gt_option_new_word("xdrop",
                              "xdrop score for extension-alignment",
                              &arguments->xdrop, (GtWord) 3);
  gt_option_parser_add_option(op, option);

  /* -opt */
  option = gt_option_new_bool("opt", "enable filtering of seeds",
                              &arguments->opt, true);
  gt_option_parser_add_option(op, option);

  /* -verbose */
  option = gt_option_new_bool("verbose", "enable verbose output",
                              &arguments->verbose, false);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_condenser_compress_arguments_check(int rest_argc,
                                                 void *tool_arguments,
                                                 GtError *err)
{
  GT_UNUSED GtCondenserCompressArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (rest_argc != 1) {
    gt_error_set(err, "give the basename of an encseq as minimum arguments");
    had_err = -1;
  }

  return had_err;
}

static int gt_condenser_compress_runner(GT_UNUSED int argc, const char **argv,
                                        int parsed_args, void *tool_arguments,
                                        GtError *err)
{
  GtCondenserCompressArguments *arguments = tool_arguments;
  GtLogger *logger;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stderr);

  if (gt_str_length(arguments->indexname) == 0UL) {
    char *basenameptr;
    basenameptr = gt_basename(argv[parsed_args]);
    gt_str_set(arguments->indexname, basenameptr);
    gt_free(basenameptr);
  }

  if (!had_err) {
    GtEncseqLoader *es_l = gt_encseq_loader_new();
    arguments->input_es = gt_encseq_loader_load(es_l, argv[parsed_args], err);
    if (arguments->input_es == NULL)
      had_err = -1;
    gt_encseq_loader_delete(es_l);
  }

  if (!had_err) {
    if (arguments->kmersize == GT_UNDEF_UINT) {
      unsigned int size =
        gt_alphabet_size(gt_encseq_alphabet(arguments->input_es));
      if (size <= 5U)
        arguments->kmersize = 10U; /* 5^10 = 9.765.625 */
      else if (size <= 7U)
        arguments->kmersize = 8U; /* 7^8  = 5.764.801 */
      else if (size <= 9U)
        arguments->kmersize = 7U; /* 8^7  = 4.782.969 */
      else if (size <= 14U)
        arguments->kmersize = 6U; /* 14^6 = 7.529.536 */
      else if (size <= 25U)
        arguments->kmersize = 5U; /* 25^5 = 9.765.625 */
      else if (size <= 55U)
        arguments->kmersize = 4U; /* 6^14 = 7.529.536 */
      else
        arguments->kmersize = 3U;
    }

    if (arguments->windowsize == GT_UNDEF_UINT) {
      arguments->windowsize = 3U * arguments->kmersize;
    }
    if (arguments->minalignlength == GT_UNDEF_UWORD) {
      arguments->minalignlength = (GtUword) (5UL * arguments->windowsize);
    }
    if (arguments->initsize == GT_UNDEF_UWORD) {
      arguments->initsize = (GtUword) (3UL * arguments->minalignlength);
    }
  }
  if (!had_err &&
      arguments->windowsize <= arguments->kmersize) {
    gt_error_set(err, "-windowsize (%u) must be larger -kmersize (%u)!",
                 arguments->windowsize, arguments->kmersize);
    had_err = -1;
  }
  if (!had_err &&
      arguments->minalignlength < (GtUword) arguments->windowsize) {
    gt_error_set(err, "-alignlength (" GT_WU ") must be at least "
                 "-windowsize (%u)!", arguments->minalignlength,
                 arguments->windowsize);
    had_err = -1;
  }
  if (!had_err && (arguments->initsize < arguments->minalignlength)) {
    gt_error_set(err, "-initsize (" GT_WU ") must be at least "
                 "-alignlength (" GT_WU ")!", arguments->initsize,
                 arguments->minalignlength);
    had_err = -1;
  }

  if (!had_err) {
    GtNREncseqCompressor *compressor;

    if (!had_err) {
      compressor = gt_n_r_encseq_compressor_new(arguments->initsize,
                                                arguments->minalignlength,
                                                arguments->max_kmer_poss,
                                                arguments->xdrop,
                                                &(arguments->scores),
                                                arguments->kmersize,
                                                arguments->windowsize,
                                                logger);
      if (!arguments->opt)
        gt_n_r_encseq_compressor_disable_opt(compressor);
      had_err = gt_n_r_encseq_compressor_compress(compressor,
                                                  arguments->indexname,
                                                  arguments->input_es, err);
      gt_n_r_encseq_compressor_delete(compressor);
    }
  }

  gt_logger_delete(logger);
  return had_err;
}

GtTool* gt_condenser_compress(void)
{
  return gt_tool_new(gt_condenser_compress_arguments_new,
                     gt_condenser_compress_arguments_delete,
                     gt_condenser_compress_option_parser_new,
                     gt_condenser_compress_arguments_check,
                     gt_condenser_compress_runner);
}
