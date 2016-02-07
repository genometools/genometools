/*
  Copyright (c) 2014 Andreas Blaufelder <9blaufel@informatik.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#include "core/alphabet_api.h"
#include "core/arraydef.h"
#include "core/basename_api.h"
#include "core/encseq_api.h"
#include "core/fa.h"
#include "core/hashmap_api.h"
#include "core/logger.h"
#include "core/logger_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/str_api.h"
#include "core/str_array_api.h"
#include "core/timer_api.h"
#include "core/unused_api.h"
#include "extended/kmer_database.h"
#include "match/sfx-mappedstr.h"
#include "tools/gt_kmer_database.h"
#include <stdio.h>
#include <unistd.h>

typedef struct {
  unsigned int kmersize;
  GtUword      sb_size,
               cutoff_value;
  bool         verbose,
               merge_only,
               cutoff,
               prune,
               mean_cutoff,
               use_hash,
               bench;
  GtStr        *print_filename;
} GtKmerDatabaseArguments;

static void* gt_kmer_database_arguments_new(void)
{
  GtKmerDatabaseArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->print_filename = gt_str_new();
  return arguments;
}

static void gt_kmer_database_arguments_delete(void *tool_arguments)
{
  GtKmerDatabaseArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->print_filename);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_kmer_database_option_parser_new(void *tool_arguments)
{
  GtKmerDatabaseArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option,
           *option_verbose,
           *option_use_cutoff,
           *option_hash,
           *option_mean_cutoff;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]",
                            "Makes a GtKmerDatabase from the input file.");

  /* -kmersize */
  option = gt_option_new_uint_min_max("kmersize", "kmersize used",
                                  &arguments->kmersize, 3U, 1U, 10U);
  gt_option_parser_add_option(op, option);

  /* -verbose */
  option_verbose = gt_option_new_bool("verbose", "prints out results of "
                                      "merging",
                                      &arguments->verbose, false);
  gt_option_parser_add_option(op, option_verbose);

  /* -merge_only */
  option = gt_option_new_bool("merge_only", "only uses merge to build DB, "
                              "doesn_t build two DBs to compare merge with a "
                              "different method (much faster). It also allows "
                              "for random intervals which are biffer than the "
                              "maximum buffer size (will be split internally).",
                              &arguments->merge_only, false);
  gt_option_parser_add_option(op, option);

  /* -use_cutoff */
  option_use_cutoff = gt_option_new_bool("use_cutoff", "uses a cutoff. see "
                                         "-set_cutoff description. Only works "
                                         "with merge_only",
                                         &arguments->cutoff, false);
  gt_option_parser_add_option(op, option_use_cutoff);
  gt_option_imply(option_use_cutoff, option);

  /* -set_cutoff */
  option = gt_option_new_uword_min("set_cutoff", "kmers occurring more often "
                                   "than this value won't be saved",
                                   &arguments->cutoff_value, (GtUword) 30,
                                   (GtUword) 1);
  gt_option_parser_add_option(op, option);
  gt_option_imply(option, option_use_cutoff);

  /* -mean_cutoff */
  option_mean_cutoff = gt_option_new_bool("mean_cutoff", "2*mean of kmer "
                                          "occurrence will be"
                                          " used as cutoff value",
                                          &arguments->mean_cutoff, false);
  gt_option_parser_add_option(op, option_mean_cutoff);
  gt_option_imply(option_mean_cutoff, option_use_cutoff);
  gt_option_exclude(option_mean_cutoff, option);

  /* -disable_prune */
  option = gt_option_new_bool("disable_prune", "disables the removel of kmers, "
                              "which occure more often than the cutoff.",
                              &arguments->prune, false);
  gt_option_parser_add_option(op, option);
  gt_option_imply(option, option_use_cutoff);

  /* -use_hash */
  option_hash = gt_option_new_bool("use_hash", "saves each kmer in kdb and "
                                   "also in a hash. afterwards both will be "
                                   "accessed and time for saving and "
                                   "accessing will be shown.",
                                   &arguments->use_hash, false);
  gt_option_parser_add_option(op, option_hash);
  gt_option_exclude(option_hash, option_use_cutoff);
  gt_option_exclude(option_hash, option_verbose);

  /* -benchmark */
  option = gt_option_new_bool("benchmark", "measures the time the tool takes to"
                              " fill the database. Doesn't test for consistency"
                              " though!", &arguments->bench, false);
  gt_option_parser_add_option(op, option);

  /* -bsize */
  option = gt_option_new_uword_min("bsize", "size of the buffer",
                                   &arguments->sb_size, (GtUword) 100000,
                                   (GtUword) 2);
  gt_option_parser_add_option(op, option);
  gt_option_exclude(option, option_hash);

  /* -outfile */
  option = gt_option_new_string("outfile", "specifies file for verbose output",
                                arguments->print_filename, NULL);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_kmer_database_arguments_check(int rest_argc,
                                            void *tool_arguments,
                                            GtError *err)
{
  GtKmerDatabaseArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (rest_argc != 1) {
    gt_error_set(err, "give the basename of an encseq");
    had_err = -1;
  }

  if (!had_err && gt_str_length(arguments->print_filename) > 0UL) {
    if (!arguments->verbose) {
      gt_error_set(err, "-outfile needs -verbose option");
      had_err = -1;
    }
  }

  return had_err;
}

GT_UNUSED static void gt_kmer_database_delete_hash_value(GtArrayGtUword *array)
{
  GT_FREEARRAY(array, GtUword);
  gt_free(array);
}

static void gt_kmer_database_add_to_hash(GtHashmap *hash, GtCodetype kmercode,
                                         GtUword position)
{
  GtArrayGtUword *arr =
    (GtArrayGtUword *) gt_hashmap_get(hash, (void *) kmercode);

  if (arr == NULL) {
    arr = gt_malloc(sizeof (*arr));
    GT_INITARRAY(arr, GtUword);
    gt_hashmap_add(hash, (void *) kmercode, (void *) arr);
  }
  if (arr->allocatedGtUword == 0)
    GT_STOREINARRAY(arr, GtUword,
                    (GtUword) 20,
                    position);
  else
    GT_STOREINARRAY(arr, GtUword,
                    arr->allocatedGtUword * 0.1,
                    position);
}

static int gt_kmer_database_runner(GT_UNUSED int argc, const char **argv,
                                   int parsed_args, void *tool_arguments,
                                   GtError *err)
{
  GtKmerDatabaseArguments *arguments = tool_arguments;
  int had_err = 0;
  GtEncseq       *es;
  GtUword        es_length,
                 nu_kmer_codes = 0;
  GtKmerDatabase *compare_db = NULL,
                 *db = NULL;
  GtLogger *logger;
  FILE *fp = NULL;
  GtHashmap *kmer_hash = NULL;
  GtTimer *timer = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  if (arguments->use_hash)
    kmer_hash = gt_hashmap_new(GT_HASH_DIRECT, NULL,
                               (GtFree) gt_kmer_database_delete_hash_value);
  if (arguments->bench)
    timer = gt_timer_new_with_progress_description("loading encoded sequence");

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stderr);

  if (arguments->verbose && gt_str_length(arguments->print_filename) > 0UL) {
    fp = gt_fa_fopen(gt_str_get(arguments->print_filename), "w", err);
    gt_logger_set_target(logger, fp);
  }

  if (!had_err) {
    GtEncseqLoader *es_l;
    if (arguments->bench)
      gt_timer_start(timer);
    es_l = gt_encseq_loader_new();
    es = gt_encseq_loader_load(es_l, argv[parsed_args], err);
    if (arguments->bench)
      gt_timer_show_progress(timer, "saving kmers (+iterating over file)",
                             stdout);
    if (es == NULL) {
      had_err = -1;
    }
    gt_encseq_loader_delete(es_l);
  }
  if (!had_err) {
    es_length = gt_encseq_total_length(es);
    if (es_length < (GtUword) arguments->kmersize) {
      gt_error_set(err, "Input is too short for used kmersize. File length: "
                   GT_WU " kmersize: %u", es_length, arguments->kmersize);
      had_err = -1;
    }
  }
  if (!had_err) {
    GtAlphabet *alphabet;
    alphabet = gt_encseq_alphabet(es);
    if (arguments->bench)
    nu_kmer_codes = gt_power_for_small_exponents(
                                            gt_alphabet_num_of_chars(alphabet),
                                            arguments->kmersize);
    if (!arguments->merge_only && !arguments->use_hash && !arguments->bench) {
      compare_db = gt_kmer_database_new(gt_alphabet_num_of_chars(alphabet),
                                        arguments->kmersize,
                                        arguments->sb_size, es);
    }
    if (!arguments->use_hash) {
      db = gt_kmer_database_new(gt_alphabet_num_of_chars(alphabet),
                                arguments->kmersize,
                                arguments->sb_size, es);
      if (arguments->cutoff) {
        if (arguments->mean_cutoff)
          gt_kmer_database_use_mean_cutoff(db, (GtUword) 2,
                                           arguments->cutoff_value);
        else
          gt_kmer_database_set_cutoff(db, arguments->cutoff_value);
        if (!arguments->prune)
          gt_kmer_database_set_prune(db);
      }
    }
  }

  if (!had_err) {
    GtUword startpos = 0,
            endpos,
            interval_id = 0;
    GtKmercodeiterator *iter;
    const GtKmercode *kmercode = NULL;
    iter = gt_kmercodeiterator_encseq_new(es, GT_READMODE_FORWARD,
                                          arguments->kmersize, 0);
    while (!had_err && startpos < es_length - (arguments->kmersize - 1)) {
      GtUword startpos_add_kmer = startpos;
      if (arguments->merge_only) {
        endpos = startpos + (arguments->kmersize - 1) +
                 (gt_rand_max((arguments->sb_size - 1) * 2));
        if (endpos > es_length)
          endpos = es_length;
      }
      else {
        endpos = startpos + (arguments->kmersize - 1) +
                 (gt_rand_max(arguments->sb_size - 1));
      }
      gt_kmercodeiterator_reset(iter, GT_READMODE_FORWARD, startpos);
      while ((kmercode = gt_kmercodeiterator_encseq_next(iter)) != NULL &&
             startpos_add_kmer <= endpos - (arguments->kmersize - 1)) {
        if (!arguments->merge_only && !arguments->use_hash &&
            !kmercode->definedspecialposition && !arguments->bench) {
          gt_kmer_database_add_kmer(compare_db, kmercode->code,
                                    startpos_add_kmer, interval_id);
        }
        if (arguments->use_hash && !kmercode->definedspecialposition) {
          gt_kmer_database_add_to_hash(kmer_hash, kmercode->code,
                                       startpos_add_kmer);
        }
        startpos_add_kmer++;
      }
      if (!arguments->use_hash) {
        gt_kmer_database_add_interval(db, startpos, endpos, interval_id++);
        gt_kmer_database_print_buffer(db, logger);
        if (!arguments->bench)
          had_err = gt_kmer_database_check_consistency(db, err);
      }
      startpos = endpos + 1;
    }
    if (!arguments->use_hash) {
      gt_kmer_database_flush(db);
      gt_kmer_database_print_buffer(db, logger);
      if (!had_err && !arguments->bench)
        had_err = gt_kmer_database_check_consistency(db, err);
      if (!arguments->merge_only && !had_err && !arguments->bench)
        had_err = gt_kmer_database_check_consistency(compare_db, err);
      if (!arguments->merge_only && !arguments->bench)
        gt_kmer_database_print(compare_db, logger, true);
      if (!arguments->merge_only && !had_err && !arguments->bench)
        had_err = gt_kmer_database_compare(compare_db, db, err);
      gt_kmer_database_print(db, logger, true);
    }
    gt_kmercodeiterator_delete(iter);
  }

  if (arguments->bench) {
    GtKmerStartpos pos;
    GtArrayGtUword *pos_hash;
    GtUword rand_access = (GtUword) 50000000,
            rand_code,
            i,
            sum = 0;
    gt_timer_show_progress(timer, "random access", stdout);
    for (i = 0; i < rand_access; i++) {
      rand_code = gt_rand_max(nu_kmer_codes - 1);
      if (arguments->use_hash) {
        pos_hash = gt_hashmap_get(kmer_hash, (const void *) rand_code);
        if (pos_hash != NULL)
          sum += pos_hash->spaceGtUword[pos_hash->nextfreeGtUword - 1];
      }
      else {
        pos = gt_kmer_database_get_startpos(db, rand_code);
        if (pos.no_positions > 0)
          sum += pos.startpos[pos.no_positions - 1];
      }
    }
    printf("sum: " GT_WU "\n", sum);

    gt_timer_show_progress(timer, "", stdout);
    gt_timer_stop(timer);
    gt_timer_delete(timer);
  }
  if (arguments->use_hash)
    gt_hashmap_delete(kmer_hash);
  gt_encseq_delete(es);
  if (!arguments->use_hash)
    gt_kmer_database_delete(db);
  if (!arguments->merge_only && !arguments->bench)
    gt_kmer_database_delete(compare_db);
  gt_logger_delete(logger);
  gt_fa_fclose(fp);

  return had_err;
}

GtTool* gt_kmer_database(void)
{
  return gt_tool_new(gt_kmer_database_arguments_new,
                     gt_kmer_database_arguments_delete,
                     gt_kmer_database_option_parser_new,
                     gt_kmer_database_arguments_check,
                     gt_kmer_database_runner);

}
