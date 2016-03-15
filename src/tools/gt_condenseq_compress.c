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
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/condenseq.h"
#include "extended/condenseq_creator.h"
#include "match/xdrop.h"
#include "tools/gt_condenseq_compress.h"

typedef struct {
  GtEncseq              *input_es;
  GtStr                 *indexname;
  GtXdropArbitraryscores scores;
  GtUword                minalignlength,
                         cutoff_value,
                         fraction,
                         initsize;
  GtWord                 xdrop;
  unsigned int           kmersize,
                         windowsize,
                         clean_percent;
  bool                   diags,
                         full_diags,
                         brute,
                         verbose,
                         kdb,
                         prune;
} GtCondenseqCompressArguments;

static void* gt_condenseq_compress_arguments_new(void)
{
  GtCondenseqCompressArguments *arguments = gt_calloc((size_t) 1,
                                                      sizeof *arguments);
  arguments->indexname = gt_str_new();
  arguments->input_es = NULL;
  return arguments;
}

static void gt_condenseq_compress_arguments_delete(void *tool_arguments)
{
  GtCondenseqCompressArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->indexname);
    gt_encseq_delete(arguments->input_es);
    gt_free(arguments);
  }
}

#define GT_CONDENSER_STR(X) #X
#define GT_CONDENSER_XSTR(X) GT_CONDENSER_STR(X)

static GtOptionParser*
gt_condenseq_compress_option_parser_new(void *tool_arguments)
{
  GtCondenseqCompressArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option,
           *option_fraction;
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
                               "length of initial unique database in bases, "
                               "should be larger than -alignlength",
                               &arguments->initsize, GT_UNDEF_UWORD);
  gt_option_parser_add_option(op, option);

  /* -alignlength */
  option = gt_option_new_uword("alignlength",
                               "required minimal length of an xdrop-alignment, "
                               "should be larger than -windowsize",
                               &arguments->minalignlength, GT_UNDEF_UWORD);
  gt_option_parser_add_option(op, option);

  /* -cutoff */
  option = gt_option_new_uword("cutoff",
                               "if a kmer is found more often than this value "
                               "it will be ignored for alignment searches. "
                               "Setting this to 0 will disable cutoffs, "
                               "leaving it undefined will use a cutoff based "
                               "on the mean number of occurrences of a k-word.",
                               &arguments->cutoff_value, GT_UNDEF_UWORD);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -fraction */
  option_fraction = gt_option_new_uword("fraction",
                               "when cutoffs aren'd disabled and no specific "
                               "value is set the mean number of occurrences "
                               "of each kmer divided by -fraction will be used "
                               "as cutoff",
                               &arguments->fraction, (GtUword) 2);
  gt_option_is_extended_option(option_fraction);
  gt_option_exclude(option, option_fraction);
  gt_option_parser_add_option(op, option_fraction);

  /* -disable_prune */
  option = gt_option_new_bool("disable_prune",
                              "when cutoffs and this option are set, "
                              "the database will still save every kmer, even "
                              "though only cutoff many kmers will be used.",
                              &arguments->prune, false);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -mat */
  option = gt_option_new_int("mat",
                             "matchscore for extension-alignment, "
                             "requirements: mat > mis, mat > 2ins, mat > 2del",
                             &arguments->scores.mat, 2);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -mis */
  option = gt_option_new_int("mis",
                             "mismatchscore for extension-alignment, ",
                             &arguments->scores.mis, -1);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -ins */
  option = gt_option_new_int("ins",
                             "insertionscore for extension-alignment",
                             &arguments->scores.ins, -2);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -del */
  option = gt_option_new_int("del",
                             "deletionscore for extension-alignment",
                             &arguments->scores.del, -2);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -xdrop */
  option = gt_option_new_word("xdrop",
                              "xdrop score for extension-alignment",
                              &arguments->xdrop, (GtWord) 3);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -brute_force */
  option = gt_option_new_bool("brute_force", "disable filtering of seeds. "
                              "Incompatible with -diagonals yes "
                              "or -full_diags yes",
                              &arguments->brute, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -diagonals */
  option = gt_option_new_bool("diagonals", "use sparse diagonals. "
                              "Incompatible with -brute_force yes. "
                              "Disabling both diagonals will result in simple "
                              "filtering of seed positions.",
                              &arguments->diags, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -full_diags */
  option = gt_option_new_bool("full_diags", "use full (time efficient "
                              "space inefficient) diagonals. "
                              "Incompatible with -brute_force yes. "
                              "Disabling both diagonals will result in simple "
                              "filtering of seed positions.",
                              &arguments->full_diags, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -clean_percent */
  option = gt_option_new_uint("diags_clean",
                              "Percentage of sparse diagonals that is allowed "
                              "to be marked as deletable. Sensible default is "
                              "set." ,
                              &arguments->clean_percent, GT_UNDEF_UINT);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -verbose */
  option = gt_option_new_bool("verbose", "enable verbose output",
                              &arguments->verbose, false);
  gt_option_parser_add_option(op, option);

  /* -kdb*/
  option = gt_option_new_bool("kdb", "prints out the kmer database (frequency "
                              "of each kmer), if -verbose each startposition "
                              "will be shown instead",
                              &arguments->kdb, false);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_condenseq_compress_arguments_check(int rest_argc,
                                                 void *tool_arguments,
                                                 GtError *err)
{
  GT_UNUSED GtCondenseqCompressArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (rest_argc != 1) {
    gt_error_set(err, "give the basename of an encseq as minimum arguments");
    had_err = -1;
  }
  if (arguments->clean_percent != GT_UNDEF_UINT &&
      arguments->clean_percent > 100U) {
    gt_error_set(err, "'-diags_clean' expects an integer argument between 0 "
                 "and 100 ");
    had_err = -1;
  }
  if (arguments->brute && (arguments->diags || arguments->full_diags)) {
    gt_error_set(err, "'-brute_force' is not compatible with '-diagonals' or "
                 "'-full_diags'");
    had_err = -1;
  }
  if (arguments->cutoff_value == 0 && arguments->prune) {
    gt_error_set(err, "'-cutoff 0' disables cutoffs, so '-disable_prune' should"
                 " not be set.");
    had_err = -1;
  }

  return had_err;
}

static int gt_condenseq_compress_runner(GT_UNUSED int argc, const char **argv,
                                        int parsed_args, void *tool_arguments,
                                        GtError *err)
{
  GtCondenseqCompressArguments *arguments = tool_arguments;
  GtLogger *logger,
           *kdb_logger;
  FILE *kmer_fp = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stderr);
  kdb_logger = gt_logger_new(arguments->kdb, GT_LOGGER_DEFLT_PREFIX, stderr);
  if (arguments->kdb) {
    kmer_fp = gt_fa_fopen("kmer_db.out", "w", err);
    gt_logger_set_target(kdb_logger, kmer_fp);
  }

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
    if (arguments->minalignlength == GT_UNDEF_UWORD)
      arguments->minalignlength = arguments->initsize != GT_UNDEF_UWORD ?
                                  arguments->initsize / (GtUword) 3UL :
                                  GT_UNDEF_UWORD;
    if (arguments->windowsize == GT_UNDEF_UINT)
      arguments->windowsize = arguments->minalignlength != GT_UNDEF_UWORD ?
                              (unsigned int) (arguments->minalignlength / 5U) :
                              GT_UNDEF_UINT;
    if (arguments->windowsize < 4U)
      arguments->windowsize = 4U;
    if (arguments->kmersize == GT_UNDEF_UINT) {
      unsigned int size =
        gt_alphabet_num_of_chars(gt_encseq_alphabet(arguments->input_es));
      /* size^k ~= 100000 */
      gt_safe_assign(arguments->kmersize,
                     gt_round_to_long(gt_log_base(100000.0, (double) size)));
      gt_logger_log(logger, "|A|: %u, k: %u",
                    size, arguments->kmersize);
    }

    if (arguments->windowsize == GT_UNDEF_UINT) {
      arguments->windowsize = 5U * arguments->kmersize;
    }
    if (arguments->minalignlength == GT_UNDEF_UWORD) {
      arguments->minalignlength = (GtUword) (3UL * arguments->windowsize);
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
    GtCondenseqCreator *ces_c;

    if (!had_err) {
      ces_c = gt_condenseq_creator_new(arguments->initsize,
                                       arguments->minalignlength,
                                       arguments->xdrop,
                                       &(arguments->scores),
                                       arguments->kmersize,
                                       arguments->windowsize,
                                       logger,
                                       err);
      if (ces_c == NULL)
        had_err = -1;
    }
    if (!had_err) {
      if (arguments->cutoff_value == GT_UNDEF_UWORD)
        gt_condenseq_creator_use_mean_cutoff(ces_c);
      else if (arguments->cutoff_value == 0)
        gt_condenseq_creator_disable_cutoff(ces_c);
      else
        gt_condenseq_creator_set_cutoff(ces_c, arguments->cutoff_value);
      gt_condenseq_creator_set_mean_fraction(ces_c, arguments->fraction);
      if (arguments->prune)
        gt_condenseq_creator_disable_prune(ces_c);
      if (arguments->brute)
        gt_condenseq_creator_enable_brute_force(ces_c);
      if (!arguments->diags)
        gt_condenseq_creator_disable_diagonals(ces_c);
      if (arguments->full_diags)
        gt_condenseq_creator_enable_full_diagonals(ces_c);
      if (arguments->clean_percent != GT_UNDEF_UINT)
        gt_condenseq_creator_set_diags_clean_limit(ces_c,
                                                   arguments->clean_percent);

      had_err = gt_condenseq_creator_create(ces_c,
                                            arguments->indexname,
                                            arguments->input_es,
                                            logger, kdb_logger, err);

      gt_condenseq_creator_delete(ces_c);
    }
  }

  gt_logger_delete(logger);
  gt_logger_delete(kdb_logger);
  if (arguments->kdb)
    gt_fa_fclose(kmer_fp);
  return had_err;
}

GtTool* gt_condenseq_compress(void)
{
  return gt_tool_new(gt_condenseq_compress_arguments_new,
                     gt_condenseq_compress_arguments_delete,
                     gt_condenseq_compress_option_parser_new,
                     gt_condenseq_compress_arguments_check,
                     gt_condenseq_compress_runner);
}
