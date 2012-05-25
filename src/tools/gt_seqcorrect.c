/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/fa.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/option_api.h"
#include "core/encseq_api.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/intbits.h"
#include "core/minmax.h"
#include "core/showtime.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread.h"
#endif
#include "tools/gt_seqcorrect.h"
#include "match/rdj-twobitenc-editor.h"
#include "match/randomcodes.h"
#include "match/randomcodes-correct.h"

#define GT_SEQCORRECT_FILESUFFIX ".cor"

typedef struct
{
  bool usefirstcodes,
       verbose,
       quiet,
       onlyaccum,
       onlyallrandomcodes,
       radixlarge,
       checksuftab,
       usemaxdepth;
  unsigned int correction_kmersize,
               samplingfactor,
               trusted_count,
               numofparts,
               radixparts,
               addbscache_depth,
               forcek;
  unsigned long maximumspace;
  GtStr *encseqinput,
        *memlimitarg;
  GtOption *refoptionmemlimit;
} GtSeqcorrectArguments;

static void* gt_seqcorrect_arguments_new(void)
{
  GtSeqcorrectArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->radixlarge = false;
  arguments->numofparts = 0;
  arguments->radixparts = 1U;
  arguments->encseqinput = gt_str_new();
  arguments->memlimitarg = gt_str_new();
  arguments->maximumspace = 0UL; /* in bytes */
  return arguments;
}

static void gt_seqcorrect_arguments_delete(void *tool_arguments)
{
  GtSeqcorrectArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->encseqinput);
  gt_option_delete(arguments->refoptionmemlimit);
  gt_str_delete(arguments->memlimitarg);
  gt_free(arguments);
}

static GtOptionParser* gt_seqcorrect_option_parser_new(void *tool_arguments)
{
  GtSeqcorrectArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionparts, *optionmemlimit, *q_option, *v_option;

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("-ii <indexname> -k <kmersize> [option ...]",
                         "K-mer based sequence correction.");

  /* -k */
  option = gt_option_new_uint_min("k", "specify the kmer size for the "
      "correction algorithm", &arguments->correction_kmersize, 31U, 2U);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -c */
  option = gt_option_new_uint_min("c", "specify the trusted count threshold",
      &arguments->trusted_count, 3U, 2U);
  gt_option_parser_add_option(op, option);

  /* -sf */
  option = gt_option_new_uint_min("sf", "specify the sampling factor",
                                  &arguments->samplingfactor, 50U, 1U);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -parts */
  optionparts = gt_option_new_uint("parts", "specify the number of parts",
                                  &arguments->numofparts, 0U);
  gt_option_parser_add_option(op, optionparts);

  /* -memlimit */
  optionmemlimit = gt_option_new_string("memlimit",
                       "specify maximal amount of memory to be used during "
                       "index construction (in bytes, the keywords 'MB' "
                       "and 'GB' are allowed)",
                       arguments->memlimitarg, NULL);
  gt_option_parser_add_option(op, optionmemlimit);
  gt_option_exclude(optionmemlimit, optionparts);
  arguments->refoptionmemlimit = gt_option_ref(optionmemlimit);

  /* -ii */
  option = gt_option_new_string("ii", "specify the input sequence",
                                arguments->encseqinput, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -forcek */
  option = gt_option_new_uint("forcek", "specify the kmersize for the bucket "
      "keys", &arguments->forcek, 0);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -usefirstcodes */
  option = gt_option_new_bool("usefirstcodes", "use first codes instead of "
      "random sampled codes", &arguments->usefirstcodes, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -v */
  v_option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, v_option);

  /* -q */
  q_option = gt_option_new_bool("q",
      "suppress standard output messages",
      &arguments->quiet, false);
  gt_option_exclude(q_option, v_option);
  gt_option_parser_add_option(op, q_option);

  /* -usemaxdepth */
  option = gt_option_new_bool("usemaxdepth", "use maxdepth in sortremaining",
                             &arguments->usemaxdepth, true);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -checksuftab */
  option = gt_option_new_bool("checksuftab", "check the suffix table",
                             &arguments->checksuftab, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -onlyaccum */
  option = gt_option_new_bool("onlyaccum", "only accumulate codes",
                             &arguments->onlyaccum, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -onlyallrandomcodes */
  option = gt_option_new_bool("onlyallrandomcodes", "only determines allcodes",
                              &arguments->onlyallrandomcodes, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -addbscachedepth */
  option = gt_option_new_uint("addbscachedepth", "only determines allcodes",
                              &arguments->addbscache_depth, 5U);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -radixlarge */
  option = gt_option_new_bool("radixlarge", "use large tables for radixsort",
                              &arguments->radixlarge, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -radixparts */
  option = gt_option_new_uint("radixparts", "specify the number of parts "
                              "for radixsort",
                              &arguments->radixparts, 1U);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  gt_option_parser_set_max_args(op, 0);

  return op;
}

static bool gt_seqcorrect_bucketkey_kmersize(GtSeqcorrectArguments *arguments,
    unsigned int *kmersize, GtError *err)
{
  bool haserr = false;
  gt_assert(kmersize != NULL);
  if (arguments->forcek > 0)
  {
    *kmersize = arguments->forcek;
    if (*kmersize >= arguments->correction_kmersize)
    {
      gt_error_set(err,"argument %u to option -forcek must be < "
          "correction kmersize (-k)", *kmersize);
      haserr = true;
    }
    else if (*kmersize > (unsigned int)GT_UNITSIN2BITENC)
    {
      gt_error_set(err,
          "argument %u to option -forcek > %u (machine word size/2)",
          *kmersize, (unsigned int)GT_UNITSIN2BITENC);
      haserr = true;
    }
  }
  else
  {
    gt_assert(arguments->correction_kmersize > 1U);
    *kmersize = MIN((unsigned int) GT_UNITSIN2BITENC,
        arguments->correction_kmersize - 1);
  }
  gt_log_log("bucketkey kmersize=%u", *kmersize);
  gt_assert(*kmersize > 0);
  return haserr;
}

static int gt_seqcorrect_arguments_check(GT_UNUSED int rest_argc,
                                         void *tool_arguments,
                                         GtError *err)
{
  GtSeqcorrectArguments *arguments = tool_arguments;
  bool haserr = false;

  gt_error_check(err);
  if (!haserr && gt_option_is_set(arguments->refoptionmemlimit))
  {
    if (gt_option_parse_spacespec(&arguments->maximumspace,
                                  "memlimit",
                                  arguments->memlimitarg,
                                  err) != 0)
    {
      haserr = true;
    }
    if (!haserr && !gt_ma_bookkeeping_enabled()) {
      gt_error_set(err, "option '-memlimit' requires "
                        "GT_MEM_BOOKKEEPING=on");
      haserr = true;
    }
  }
#ifdef GT_THREADS_ENABLED
  if (!haserr) {
    if (gt_jobs > 1 && gt_ma_bookkeeping_enabled()) {
      gt_error_set(err, "gt option '-j' and GT_MEM_BOOKKEEPING=on "
                        "are incompatible");
      haserr = true;
    }
  }
#endif
  return haserr ? -1 : 0;
}

static int gt_seqcorrect_apply_corrections(GtEncseq *encseq,
    const char *indexname, const unsigned int threads,
    GtError *err)
{
  bool haserr = false;
  GtTwobitencEditor *editor;
  unsigned int threadcount;
  editor = gt_twobitenc_editor_new(encseq, indexname, err);
  if (editor == NULL)
    haserr = true;
  for (threadcount = 0; !haserr && threadcount < threads; threadcount++)
  {
    FILE *corrections;
    GtStr *filename;
    filename = gt_str_new_cstr(indexname);
    gt_str_append_char(filename, '.');
    gt_str_append_uint(filename, threadcount);
    gt_str_append_cstr(filename, GT_SEQCORRECT_FILESUFFIX);
    corrections = gt_fa_fopen(gt_str_get(filename), "r", err);
    if (corrections == NULL)
      haserr = true;
    else
    {
      unsigned long pos;
      GtUchar newchar;
      size_t retval;
      while ((retval = fread(&pos, sizeof (pos), (size_t)1, corrections))
          == (size_t)1)
      {
        newchar = (GtUchar)(pos & 3UL);
        pos >>= 2;
        gt_twobitenc_editor_edit(editor, pos, newchar);
      }
      if (ferror(corrections) != 0)
      {
        gt_error_set(err, "error by reading file %s", gt_str_get(filename));
        haserr = true;
      }
      gt_fa_fclose(corrections);
    }
    gt_str_delete(filename);
  }
  gt_twobitenc_editor_delete(editor);
  return haserr ? -1 : 0;
}

static int gt_seqcorrect_runner(GT_UNUSED int argc,
                                GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments,
                                GtError *err)
{
  GtSeqcorrectArguments *arguments = tool_arguments;
  GtEncseqLoader *el = NULL;
  GtEncseq *encseq = NULL;
  GtTimer *timer = NULL;
  GtLogger *default_logger, *verbose_logger;
  bool haserr = false;

  gt_error_check(err);
  gt_assert(arguments);

  if (gt_showtime_enabled())
  {
    timer = gt_timer_new_with_progress_description("for initialization");
    gt_timer_show_cpu_time_by_progress(timer);
    gt_timer_start(timer);
  }
  default_logger =
    gt_logger_new(!arguments->quiet, GT_LOGGER_DEFLT_PREFIX, stdout);
  gt_logger_log(default_logger, "gt seqcorrect");
  verbose_logger =
    gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);

  el = gt_encseq_loader_new();
  gt_encseq_loader_drop_description_support(el);
  gt_encseq_loader_disable_autosupport(el);
  encseq = gt_encseq_loader_load(el, gt_str_get(arguments->encseqinput),
      err);
  if (encseq == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (gt_encseq_mirror(encseq, err) != 0)
    {
      haserr = true;
    }
  }

  if (!haserr)
  {
    GtRandomcodesCorrectData **data_array = NULL;
    unsigned int bucketkey_kmersize, threadcount;
    unsigned long nofkmergroups = 0, nofkmeritvs = 0, nofcorrections = 0,
                  nofkmers = 0;

#ifdef GT_THREADS_ENABLED
    const unsigned int threads = gt_jobs;
#else
    const unsigned int threads = 1U;
#endif

    data_array = gt_malloc(sizeof (*data_array) * threads);
    for (threadcount = 0; !haserr && threadcount < threads; threadcount++)
    {
      data_array[threadcount] = gt_randomcodes_correct_data_new(encseq,
          arguments->correction_kmersize, arguments->trusted_count,
          gt_str_get(arguments->encseqinput), GT_SEQCORRECT_FILESUFFIX,
          threadcount, err);
      if ((data_array[threadcount]) == NULL)
      {
        haserr = true;
      }
    }
    gt_log_log("correction kmersize=%u", arguments->correction_kmersize);
    haserr = gt_seqcorrect_bucketkey_kmersize(arguments,
        &bucketkey_kmersize, err);
    if (!haserr)
    {
      if (storerandomcodes_getencseqkmers_twobitencoding(encseq,
            bucketkey_kmersize,
            arguments->numofparts,
            arguments->maximumspace,
            arguments->correction_kmersize,
            arguments->usefirstcodes,
            arguments->samplingfactor,
            arguments->usemaxdepth,
            arguments->checksuftab,
            arguments->onlyaccum,
            arguments->onlyallrandomcodes,
            arguments->addbscache_depth,
            0,
            arguments->radixlarge ? false : true,
            arguments->radixparts,
            gt_randomcodes_correct_process_bucket,
            NULL,
            data_array,
            verbose_logger,
            timer,
            err) != 0)
            {
              haserr = true;
            }
    }
    for (threadcount = 0; threadcount < threads; threadcount++)
    {
      if (!haserr)
        gt_randomcodes_correct_data_collect_stats(data_array[threadcount],
            threadcount, &nofkmergroups, &nofkmeritvs, &nofkmers,
            &nofcorrections);
      gt_randomcodes_correct_data_delete(data_array[threadcount]);
    }
    gt_logger_log(verbose_logger, "total number of k-mers: %lu",
        nofkmers);
    gt_logger_log(verbose_logger, "number of different k-mers: %lu",
        nofkmeritvs);
    gt_logger_log(verbose_logger, "number of different k-1-mers: %lu",
        nofkmergroups);
    gt_logger_log(verbose_logger, "number of kmer corrections: %lu",
        nofcorrections);
    gt_free(data_array);
    if (!haserr)
      gt_seqcorrect_apply_corrections(encseq,
          gt_str_get(arguments->encseqinput), threads, err);
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  if (timer != NULL)
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
  gt_logger_delete(default_logger);
  gt_logger_delete(verbose_logger);
  return haserr ? -1 : 0;
}

GtTool* gt_seqcorrect(void)
{
  return gt_tool_new(gt_seqcorrect_arguments_new,
                     gt_seqcorrect_arguments_delete,
                     gt_seqcorrect_option_parser_new,
                     gt_seqcorrect_arguments_check,
                     gt_seqcorrect_runner);
}
