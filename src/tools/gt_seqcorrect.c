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
#include "core/undef_api.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif
#include "core/warning_api.h"
#include "core/xposix.h"
#include "tools/gt_seqcorrect.h"
#include "match/reads2twobit.h"
#include "match/rdj-contfinder.h"
#include "match/rdj-twobitenc-editor.h"
#include "match/randomcodes.h"
#include "match/randomcodes-correct.h"
#include "tools/gt_seqcorrect.h"

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
        *memlimitarg,
        *indexname;
  GtOption *refoptionmemlimit;
  GtStrArray *db;
  bool phred64;
  unsigned long maxlow;
  unsigned int lowqual;
} GtSeqcorrectArguments;

static void* gt_seqcorrect_arguments_new(void)
{
  GtSeqcorrectArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->radixlarge = false;
  arguments->numofparts = 0;
  arguments->radixparts = 1U;
  arguments->encseqinput = gt_str_new();
  arguments->indexname = gt_str_new();
  arguments->memlimitarg = gt_str_new();
  arguments->maximumspace = 0UL; /* in bytes */
  arguments->db = gt_str_array_new();
  return arguments;
}

static void gt_seqcorrect_arguments_delete(void *tool_arguments)
{
  GtSeqcorrectArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->encseqinput);
  gt_str_delete(arguments->indexname);
  gt_option_delete(arguments->refoptionmemlimit);
  gt_str_delete(arguments->memlimitarg);
  gt_str_array_delete(arguments->db);
  gt_free(arguments);
}

static GtOptionParser* gt_seqcorrect_option_parser_new(void *tool_arguments)
{
  GtSeqcorrectArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionparts, *optionmemlimit, *q_option, *v_option,
           *db_option, *maxlow_option;

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("(-ii <indexname>|-db <filenames>) "
      "-k <kmersize> [option ...]", "K-mer based sequence correction.");

  /* -db */
  db_option = gt_option_new_filename_array("db",
      GT_READS2TWOBIT_LIBSPEC_HELPMSG, arguments->db);
  gt_option_hide_default(db_option);
  gt_option_parser_add_option(op, db_option);

  /* -indexname */
  option = gt_option_new_string("indexname", "specify the indexname to use "
      "for the input\ndefault: first argument of the -db option",
      arguments->indexname, NULL);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* -ii */
  option = gt_option_new_string("ii", "specify the input enseq index",
                                arguments->encseqinput, NULL);
  gt_option_is_mandatory_either(option, db_option);
  gt_option_exclude(option, db_option);
  gt_option_parser_add_option(op, option);

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

  /* -maxlow */
  maxlow_option = gt_option_new_ulong("maxlow",
      "maximal number of low-quality positions in a read\n"
      "default: infinite",
      &arguments->maxlow, GT_UNDEF_ULONG);
  gt_option_hide_default(maxlow_option);
  gt_option_is_extended_option(maxlow_option);
  gt_option_parser_add_option(op, maxlow_option);

  /* -lowqual */
  option = gt_option_new_uint_max("lowqual",
      "maximal quality for a position to be considered low-quality",
      &arguments->lowqual, 3U, 127U);
  gt_option_is_extended_option(option);
  gt_option_imply(option, maxlow_option);
  gt_option_parser_add_option(op, option);

  /* -phred64 */
  option = gt_option_new_bool("phred64", "use phred64 scores for FastQ format",
      &arguments->phred64, false);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

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
  gt_log_log("number of correction lists: %u", threads);
  for (threadcount = 0; !haserr && threadcount < threads; threadcount++)
  {
    FILE *corrections;
    GtStr *filename;
    filename = gt_str_new_cstr(indexname);
    gt_str_append_char(filename, '.');
    gt_str_append_uint(filename, threadcount);
    gt_str_append_cstr(filename, GT_SEQCORRECT_FILESUFFIX);
    corrections = gt_fa_fopen(gt_str_get(filename), "r", err);
    gt_log_log("apply corrections list %s.%u.%s", indexname, threadcount,
        GT_SEQCORRECT_FILESUFFIX);
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
  gt_logger_log(verbose_logger, "verbose output enabled");

  if (gt_str_array_size(arguments->db) != 0)
  {
    GtReads2Twobit *r2t;
    unsigned long i;
    bool autoindexname = false;
    gt_logger_log(verbose_logger, "input is a list of libraries");
    if (gt_str_length(arguments->indexname) == 0)
    {
      autoindexname = true;
      gt_str_append_cstr(arguments->indexname,
          gt_str_array_get(arguments->db, 0));
    }
    gt_logger_log(verbose_logger, "indexname: %s [%s]",
          gt_str_get(arguments->indexname), autoindexname ?
          "first input library" : "user-specified");

    r2t = gt_reads2twobit_new(arguments->indexname);

    if (arguments->phred64)
      gt_reads2twobit_use_phred64(r2t);

    if (arguments->maxlow != GT_UNDEF_ULONG)
      gt_reads2twobit_set_quality_filter(r2t, arguments->maxlow,
          (char)arguments->lowqual);

    for (i = 0; i < gt_str_array_size(arguments->db) && !haserr; i++)
    {
      GtStr *dbentry = gt_str_array_get_str(arguments->db, i);
      if (gt_reads2twobit_add_library(r2t, dbentry, err) != 0)
        haserr = true;
    }
    if (!haserr)
    {
      if (gt_reads2twobit_encode(r2t, err) != 0)
        haserr = true;
    }
    gt_assert(gt_str_length(arguments->encseqinput) == 0);
    gt_str_append_cstr(arguments->encseqinput,
        gt_str_get(arguments->indexname));
    if (!haserr)
    {
      unsigned long nofreads_valid, nofreads_invalid, nofreads_input,
                    tlen_valid;
      bool varlen;
      nofreads_valid = gt_reads2twobit_nofseqs(r2t);
      nofreads_invalid = gt_reads2twobit_nof_invalid_seqs(r2t);
      nofreads_input = nofreads_valid + nofreads_invalid;
      tlen_valid = gt_reads2twobit_total_seqlength(r2t) -
        gt_reads2twobit_nofseqs(r2t);

      gt_logger_log(default_logger, "number of reads in original read set "
          "= %lu", nofreads_input);

      varlen = (gt_reads2twobit_seqlen_eqlen(r2t) == 0);
      if (varlen)
        gt_logger_log(verbose_logger, "read length = variable [%lu..%lu]",
            gt_reads2twobit_seqlen_min(r2t), gt_reads2twobit_seqlen_max(r2t));
      else
        gt_logger_log(verbose_logger, "read length = %lu",
            gt_reads2twobit_seqlen_eqlen(r2t) - 1UL);

      gt_logger_log(verbose_logger, "total length of original read set = %lu",
          tlen_valid + gt_reads2twobit_invalid_seqs_totallength(r2t));
      gt_logger_log(verbose_logger, "low-quality reads = %lu "
          "[%.2f %% of input]", nofreads_invalid, (float)nofreads_invalid *
          100 / (float)nofreads_input);
      if (!arguments->verbose)
        gt_logger_log(default_logger, "low-quality reads = %lu",
            nofreads_invalid);
      if (!haserr)
      {
        if (gt_reads2twobit_write_encseq(r2t, err) != 0)
          haserr = true;
        if (!haserr)
        {
          gt_logger_log(verbose_logger,
              "number of reads in output read set = %lu", nofreads_valid);
          gt_logger_log(verbose_logger,
              "total length of output read set = %lu", tlen_valid);
          gt_logger_log(verbose_logger, "read set saved as GtEncseq: %s.%s",
              gt_str_get(arguments->indexname), varlen ?
              "(esq|ssp)" : "esq");
        }
      }
    }
    gt_reads2twobit_delete(r2t);
  }
  else
  {
    gt_logger_log(verbose_logger, "input is an encseq: %s",
        gt_str_get(arguments->encseqinput));
  }
  if (!haserr)
  {
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
      if (!haserr && gt_encseq_has_md5_support(encseq))
      {
        GtStr *md5path = gt_str_clone(arguments->encseqinput);
        gt_str_append_cstr(md5path, GT_MD5TABFILESUFFIX);
        gt_xunlink(gt_str_get(md5path));
        gt_warning("MD5 support detected -- sequence correction will "
                   "invalidate MD5 table, permanently disabling MD5 support "
                   "in index %s", gt_str_get(arguments->encseqinput));
        gt_str_delete(md5path);
      }
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
      if (!haserr) {
        gt_randomcodes_correct_data_collect_stats(data_array[threadcount],
            threadcount, &nofkmergroups, &nofkmeritvs, &nofkmers,
            &nofcorrections);
      }
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
    if (!haserr) {
      if (gt_seqcorrect_apply_corrections(encseq,
          gt_str_get(arguments->encseqinput), threads, err) != 0) {
        haserr = true;
      }
    }
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
