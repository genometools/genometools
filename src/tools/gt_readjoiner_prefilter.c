/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/log.h"
#include "core/logger.h"
#include "core/encseq.h"
#include "core/fa.h"
#include "core/ma.h"
#include "core/output_file_api.h"
#include "core/option_api.h"
#include "core/unused_api.h"
#include "core/undef_api.h"
#include "match/rdj-contfinder.h"
#include "match/rdj-filesuf-def.h"
#include "match/rdj-version.h"
#include "match/reads_library.h"
#include "match/reads2twobit.h"
#include "tools/gt_readjoiner_prefilter.h"

typedef struct {
  bool verbose, quiet;
  bool singlestrand, encodeonly, cntlist, encseq, seqnums, fasta, copynum,
       libtable, phred64;
  GtStr *readset;
  GtStrArray *db;
  /* rdj-radixsort test */
  bool testrs, testrs_print;
  unsigned long testrs_depth, testrs_maxdepth, maxlow;
  unsigned int lowqual;
} GtReadjoinerPrefilterArguments;

static void* gt_readjoiner_prefilter_arguments_new(void)
{
  GtReadjoinerPrefilterArguments *arguments = gt_calloc((size_t)1,
      sizeof *arguments);
  arguments->readset = gt_str_new();
  arguments->db = gt_str_array_new();
  return arguments;
}

static void gt_readjoiner_prefilter_arguments_delete(void *tool_arguments)
{
  GtReadjoinerPrefilterArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->readset);
  gt_str_array_delete(arguments->db);
  gt_free(arguments);
}

static GtOptionParser* gt_readjoiner_prefilter_option_parser_new(
    void *tool_arguments)
{
  GtReadjoinerPrefilterArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *singlestrand_option, *encodeonly_option, *cntlist_option,
           *fasta_option, *seqnums_option, *encseq_option, *readset_option,
           *v_option, *q_option, *db_option, *copynum_option, *libtable_option,
           *testrs_option, *testrs_depth_option, *testrs_print_option,
           *testrs_maxdepth_option, *phred64_option, *maxlow_option,
           *lowqual_option;

  gt_assert(arguments);

  op = gt_option_parser_new("[option ...]",
                            "Remove contained and low-quality reads "
                            "and encode read set in GtEncseq format.");

  /* -readset */
  readset_option = gt_option_new_string("readset",
      "specify the readset name\n"
      "default: filename of first input sequence_file",
      arguments->readset, NULL);
  gt_option_hide_default(readset_option);
  gt_option_parser_add_option(op, readset_option);

  /* -db */
  db_option = gt_option_new_filename_array("db",
      GT_READS2TWOBIT_LIBSPEC_HELPMSG, arguments->db);
  gt_option_hide_default(db_option);
  gt_option_is_mandatory(db_option);
  gt_option_parser_add_option(op, db_option);

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
  lowqual_option = gt_option_new_uint_max("lowqual",
      "maximal quality for a position to be considered low-quality",
      &arguments->lowqual, 3U, 127U);
  gt_option_is_extended_option(lowqual_option);
  gt_option_imply(lowqual_option, maxlow_option);
  gt_option_parser_add_option(op, lowqual_option);

  /* -phred64 */
  phred64_option = gt_option_new_bool("phred64",
      "use phred64 scores for FastQ format",
      &arguments->phred64, false);
  gt_option_is_extended_option(phred64_option);
  gt_option_parser_add_option(op, phred64_option);

  /* -singlestrand */
  singlestrand_option = gt_option_new_bool("singlestrand",
      "do not use reverse complements of the reads",
      &arguments->singlestrand, false);
  gt_option_is_development_option(singlestrand_option);
  gt_option_parser_add_option(op, singlestrand_option);

  /* -encodeonly */
  encodeonly_option = gt_option_new_bool("encodeonly",
      "do not remove contained reads",
      &arguments->encodeonly, false);
  gt_option_is_development_option(encodeonly_option);
  gt_option_parser_add_option(op, encodeonly_option);

  /* -encseq */
  encseq_option = gt_option_new_bool("encseq",
      "output reads in GtEncseq format",
      &arguments->encseq, true);
  gt_option_is_development_option(encseq_option);
  gt_option_parser_add_option(op, encseq_option);

  /* -libtable */
  libtable_option = gt_option_new_bool("libtable",
      "output reads libraries table",
      &arguments->libtable, true);
  gt_option_is_development_option(libtable_option);
  gt_option_parser_add_option(op, libtable_option);

  /* -fasta */
  fasta_option = gt_option_new_bool("fasta",
      "output reads in MultiFasta format",
      &arguments->fasta, false);
  gt_option_is_development_option(fasta_option);
  gt_option_parser_add_option(op, fasta_option);

  /* -cnt */
  cntlist_option = gt_option_new_bool("cnt",
      "output contained reads list",
      &arguments->cntlist, false);
  gt_option_is_development_option(cntlist_option);
  gt_option_exclude(encodeonly_option, cntlist_option);
  gt_option_parser_add_option(op, cntlist_option);

  /* -seqnums */
  seqnums_option = gt_option_new_bool("seqnums",
      "output sorted sequence numbers",
      &arguments->seqnums, false);
  gt_option_is_development_option(seqnums_option);
  gt_option_exclude(encodeonly_option, seqnums_option);
  gt_option_parser_add_option(op, seqnums_option);

  /* -copynum */
  copynum_option = gt_option_new_bool("copynum",
      "[eqlen only] output reads copy number to <readset>"
      GT_READJOINER_SUFFIX_READSCOPYNUM,
      &arguments->copynum, false);
  gt_option_is_development_option(copynum_option);
  gt_option_exclude(encodeonly_option, copynum_option);
  gt_option_parser_add_option(op, copynum_option);

  /* -testrs */
  testrs_option = gt_option_new_bool("testrs",
      "run gt_radixsort_str test (match/radixsort_str.[ch])",
      &arguments->testrs, false);
  gt_option_is_development_option(testrs_option);
  gt_option_parser_add_option(op, testrs_option);

  /* -testrs-print */
  testrs_print_option = gt_option_new_bool("testrs-print",
      "printf gt_radixsort_str test results",
      &arguments->testrs_print, true);
  gt_option_is_development_option(testrs_print_option);
  gt_option_parser_add_option(op, testrs_print_option);

  /* -testrs-depth */
  testrs_depth_option = gt_option_new_ulong("testrs-depth",
      "depth for gt_radixsort_str test",
      &arguments->testrs_depth, 0);
  gt_option_is_development_option(testrs_depth_option);
  gt_option_parser_add_option(op, testrs_depth_option);

  /* -testrs-maxdepth */
  testrs_maxdepth_option = gt_option_new_ulong("testrs-maxdepth",
      "depth for gt_radixsort_str test",
      &arguments->testrs_maxdepth, 0);
  gt_option_is_development_option(testrs_maxdepth_option);
  gt_option_parser_add_option(op, testrs_maxdepth_option);

  gt_option_parser_set_version_func(op, gt_readjoiner_show_version);
  gt_option_parser_set_max_args(op, 0U);
  return op;
}

static void gt_readjoiner_prefilter_list_input_files(
    GtReadjoinerPrefilterArguments *arguments, GtLogger *verbose_logger)
{
  unsigned long i;
  GtStr *inputfileslist = gt_str_new();
  for (i = 0; i < gt_str_array_size(arguments->db); i++)
  {
    gt_str_append_cstr(inputfileslist, gt_str_array_get(arguments->db, i));
    if (i + 1 < gt_str_array_size(arguments->db))
      gt_str_append_cstr(inputfileslist, ", ");
  }
  gt_logger_log(verbose_logger, "input files = %s", gt_str_get(inputfileslist));
  gt_str_delete(inputfileslist);
}

#define GT_READJOINER_PREFILTER_CF_OUTPUT(FUNC, FILESUF, MSG)\
{\
  GtStr *fn;\
  fn = gt_str_new_cstr(gt_str_get(arguments->readset));\
  gt_str_append_cstr(fn, (FILESUF));\
  had_err = (FUNC)(contfinder, gt_str_get(fn), err);\
  gt_logger_log(verbose_logger, MSG": %s", gt_str_get(fn));\
  gt_str_delete(fn);\
}

static int gt_readjoiner_prefilter_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    void *tool_arguments, GtError *err)
{
  GtReadjoinerPrefilterArguments *arguments = tool_arguments;
  int had_err = 0;
  bool varlen;
  unsigned long i;
  unsigned long nofreads_valid, nofreads_invalid, nofreads_input,
                nofreads_output,
                tlen_valid, tlen_invalid, tlen_input;
  GtLogger *default_logger, *verbose_logger;
  GtReads2Twobit *r2t;

  default_logger =
    gt_logger_new(!arguments->quiet, GT_LOGGER_DEFLT_PREFIX, stdout);
  verbose_logger =
    gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);

  gt_logger_log(default_logger,
      "gt readjoiner prefilter (version "GT_READJOINER_VERSION")");

  /* default readset name: first db argument */
  if (gt_str_length(arguments->readset) == 0)
    gt_str_append_cstr(arguments->readset, gt_str_array_get(arguments->db, 0));
  gt_logger_log(verbose_logger, "readset name = %s",
      gt_str_get(arguments->readset));

  if (arguments->verbose)
    gt_readjoiner_prefilter_list_input_files(arguments, verbose_logger);

  r2t = gt_reads2twobit_new(arguments->readset);

  if (arguments->phred64)
    gt_reads2twobit_use_phred64(r2t);

  if (arguments->maxlow != GT_UNDEF_ULONG)
    gt_reads2twobit_set_quality_filter(r2t, arguments->maxlow,
        (char)arguments->lowqual);

  for (i = 0; i < gt_str_array_size(arguments->db) && !had_err; i++)
  {
    GtStr *dbentry = gt_str_array_get_str(arguments->db, i);
    had_err = gt_reads2twobit_add_library(r2t, dbentry, err);
  }

  if (!had_err)
    had_err = gt_reads2twobit_encode(r2t, err);

  if (!had_err)
  {
    nofreads_valid = gt_reads2twobit_nofseqs(r2t);
    nofreads_invalid = gt_reads2twobit_nof_invalid_seqs(r2t);
    nofreads_input = nofreads_valid + nofreads_invalid;
    tlen_valid = gt_reads2twobit_total_seqlength(r2t) -
      gt_reads2twobit_nofseqs(r2t);
    tlen_invalid = gt_reads2twobit_invalid_seqs_totallength(r2t);
    tlen_input = tlen_valid + tlen_invalid;

    gt_logger_log(default_logger, "number of reads in complete readset = %lu",
        nofreads_input);

    varlen = (gt_reads2twobit_seqlen_eqlen(r2t) == 0);
    if (varlen)
      gt_logger_log(verbose_logger, "read length = variable [%lu..%lu]",
          gt_reads2twobit_seqlen_min(r2t), gt_reads2twobit_seqlen_max(r2t));
    else
      gt_logger_log(verbose_logger, "read length = %lu",
          gt_reads2twobit_seqlen_eqlen(r2t) - 1UL);

    gt_logger_log(verbose_logger, "total length of complete readset = %lu",
        tlen_input);
    gt_logger_log(verbose_logger, "low-quality reads = %lu "
        "[%.2f %% of input]", nofreads_invalid, (float)nofreads_invalid * 100 /
        (float)nofreads_input);
    if (!arguments->verbose)
      gt_logger_log(default_logger, "low-quality reads = %lu",
          nofreads_invalid);
    nofreads_output = nofreads_valid;
    if (arguments->encodeonly)
    {
      if (!had_err && arguments->fasta)
      {
        GtStr *fn;
        fn = gt_str_new_cstr(gt_str_get(arguments->readset));
        gt_str_append_cstr(fn, GT_READJOINER_SUFFIX_PREFILTERED_FAS);
        had_err = gt_reads2twobit_write_fasta(r2t, gt_str_get(fn), NULL, err);
        gt_logger_log(verbose_logger, "readset saved: %s", gt_str_get(fn));
        gt_str_delete(fn);
      }
      if (!had_err && arguments->encseq)
      {
        had_err = gt_reads2twobit_write_encseq(r2t, err);
        gt_logger_log(verbose_logger, "readset saved: %s.%s",
            gt_str_get(arguments->readset), varlen ? "(esq|ssp)" : "esq");
      }
      gt_logger_log(default_logger, "number of reads in output readset = %lu",
          nofreads_output);
    }
    else
    {
      GtContfinder *contfinder;
      unsigned long nofreads_contained;
      contfinder = gt_contfinder_new(r2t);
      gt_contfinder_run(contfinder, !arguments->singlestrand,
          arguments->copynum);

      nofreads_contained = gt_contfinder_nofcontained(contfinder);
      nofreads_output -= nofreads_contained;

      gt_logger_log(verbose_logger, "contained reads = %lu [%.2f %% of input]",
          nofreads_contained, (float)nofreads_contained * 100 /
          (float)nofreads_input);
      if (!arguments->verbose)
        gt_logger_log(default_logger, "contained reads = %lu",
            nofreads_contained);

      gt_logger_log(default_logger, "number of reads in filtered readset = %lu",
          nofreads_output);

      if (!had_err && arguments->cntlist)
      {
        GT_READJOINER_PREFILTER_CF_OUTPUT(gt_contfinder_write_cntlist,
            GT_READJOINER_SUFFIX_CNTLIST, "contained reads list saved");
      }
      if (!had_err && arguments->copynum)
      {
        GT_READJOINER_PREFILTER_CF_OUTPUT(gt_contfinder_write_copynum,
            GT_READJOINER_SUFFIX_READSCOPYNUM, "reads copy number saved");
      }
      if (!had_err && arguments->seqnums)
      {
        GT_READJOINER_PREFILTER_CF_OUTPUT(gt_contfinder_write_sorted_seqnums,
            GT_READJOINER_SUFFIX_SEQNUMS, "sorted sequence numbers "
            "saved");
      }
      if (!had_err && arguments->fasta)
      {
        GtStr *fn;
        fn = gt_str_new_cstr(gt_str_get(arguments->readset));
        gt_str_append_cstr(fn, GT_READJOINER_SUFFIX_PREFILTERED_FAS);
        had_err = gt_reads2twobit_write_fasta(r2t, gt_str_get(fn),
            gt_contfinder_contained(contfinder), err);
        gt_logger_log(verbose_logger, "suffix-prefix-free readset saved: %s",
            gt_str_get(fn));
        gt_str_delete(fn);
      }
      if (!had_err && arguments->encseq)
      {
        gt_reads2twobit_delete_sequences(r2t,
            gt_contfinder_contained(contfinder));
        had_err = gt_reads2twobit_write_encseq(r2t, err);
        gt_logger_log(verbose_logger, "suffix-prefix-free readset saved: %s.%s",
            gt_str_get(arguments->readset),
            varlen ? "(esq|ssp)" : "esq");
      }
      if (arguments->testrs)
      {
        gt_contfinder_radixsort_str_eqlen_tester(contfinder,
            !arguments->singlestrand,
            arguments->testrs_depth, arguments->testrs_maxdepth,
            arguments->testrs_print);
      }
      gt_contfinder_delete(contfinder);
    }
    if (!had_err && arguments->libtable) {
      GtStr *fn;
      fn = gt_str_new_cstr(gt_str_get(arguments->readset));
      gt_str_append_cstr(fn, GT_READS_LIBRARY_TABLE_FILESUFFIX);
      had_err = gt_reads2twobit_write_libraries_table(r2t, gt_str_get(fn), err);
      gt_logger_log(verbose_logger, "reads library table saved: %s",
          gt_str_get(fn));
      gt_str_delete(fn);
    }
  }
  gt_reads2twobit_delete(r2t);
  gt_logger_delete(verbose_logger);
  gt_logger_delete(default_logger);
  return had_err;
}

GtTool* gt_readjoiner_prefilter(void)
{
  return gt_tool_new(gt_readjoiner_prefilter_arguments_new,
                     gt_readjoiner_prefilter_arguments_delete,
                     gt_readjoiner_prefilter_option_parser_new,
                     NULL,
                     gt_readjoiner_prefilter_runner);
}
