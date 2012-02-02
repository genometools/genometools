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
#include "core/outputfile.h"
#include "core/option_api.h"
#include "core/unused_api.h"
#include "match/rdj-contfinder.h"
#include "match/rdj-filesuf-def.h"
#include "tools/gt_readjoiner_prefilter.h"

typedef struct {
  bool verbose, quiet;
  bool singlestrand, encodeonly, cntlist, encseqall, encseq, sorted, seqnums,
       fasta, twobit, seppos, testrs;
  unsigned long width;
  GtStr *readset;
  GtStrArray *db;
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
           *fasta_option, *sorted_option, *seqnums_option, *twobit_option,
           *seppos_option, *testrs_option, *encseqall_option,
           *encseq_option, *readset_option, *v_option, *q_option, *db_option;

  gt_assert(arguments);

  op = gt_option_parser_new("[option ...]",
                            "Remove contained reads and reads containing "
                            "ambiguity codes in given sequence_file(s).");

  /* -readset */
  readset_option = gt_option_new_string("readset",
      "specify the readset name\n"
      "default: filename of first input sequence_file",
      arguments->readset, NULL);
  gt_option_hide_default(readset_option);
  gt_option_parser_add_option(op, readset_option);

  /* -db */
  db_option = gt_option_new_filename_array("db","specify the name of the "
      "input files (Fasta format)", arguments->db);
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

  /* -singlestrand */
  singlestrand_option = gt_option_new_bool("singlestrand",
      "do not use reverse complements of the reads",
      &arguments->singlestrand, false);
  gt_option_is_development_option(singlestrand_option);
  gt_option_parser_add_option(op, singlestrand_option);

  /* -sorted */
  sorted_option = gt_option_new_bool("sorted",
      "output sorted sequences/sequence numbers\n"
      "contained reads are not eliminated",
      &arguments->sorted, false);
  gt_option_is_development_option(sorted_option);
  gt_option_parser_add_option(op, sorted_option);

  /* -encseq */
  encseq_option = gt_option_new_bool("encseq",
      "output in encseq format",
      &arguments->encseq, true);
  gt_option_is_development_option(encseq_option);
  gt_option_parser_add_option(op, encseq_option);

  /* -fasta */
  fasta_option = gt_option_new_bool("fasta",
      "output in fasta format",
      &arguments->fasta, false);
  gt_option_is_development_option(fasta_option);
  gt_option_parser_add_option(op, fasta_option);

  /* -seqnums */
  seqnums_option = gt_option_new_bool("seqnums",
      "output sequence numbers",
      &arguments->seqnums, false);
  gt_option_is_development_option(seqnums_option);
  gt_option_parser_add_option(op, seqnums_option);

  /* -2bit */
  twobit_option = gt_option_new_bool("2bit",
      "output 2bit representation",
      &arguments->twobit, false);
  gt_option_is_development_option(twobit_option);
  gt_option_parser_add_option(op, twobit_option);

  /* -cnt */
  cntlist_option = gt_option_new_bool("cnt",
      "output contained reads list",
      &arguments->cntlist, false);
  gt_option_is_development_option(cntlist_option);
  gt_option_parser_add_option(op, cntlist_option);

  /* -seppos */
  seppos_option = gt_option_new_bool("seppos",
      "output separator positions",
      &arguments->seppos, false);
  gt_option_is_development_option(seppos_option);
  gt_option_parser_add_option(op, seppos_option);

  /* -encodeonly */
  encodeonly_option = gt_option_new_bool("encodeonly",
      "exit after parsing/encoding step (before sorting)",
      &arguments->encodeonly, false);
  gt_option_is_development_option(encodeonly_option);
  gt_option_parser_add_option(op, encodeonly_option);

  /* -testrs */
  testrs_option = gt_option_new_bool("testrs",
      "run gt_radixsort test (match/radixsort.[ch])",
      &arguments->testrs, false);
  gt_option_is_development_option(testrs_option);
  gt_option_parser_add_option(op, testrs_option);

  /* -encseqall */
  encseqall_option = gt_option_new_bool("encseqall",
      "output encseq before containments removal",
      &arguments->encseqall, false);
  gt_option_is_development_option(encseqall_option);
  gt_option_parser_add_option(op, encseqall_option);

  gt_option_parser_set_max_args(op, 0U);
  return op;
}

static int gt_prefilter_output_varlen_encseq(const char *indexname,
    GtError *err)
{
  GtStr *fas_path;
  GtEncseqEncoder *encoder = gt_encseq_encoder_new();
  GtStrArray *infiles = gt_str_array_new();
  int had_err = 0;

  fas_path = gt_str_new_cstr(indexname);
  gt_str_append_cstr(fas_path, GT_READJOINER_SUFFIX_PREFILTERED_FAS);
  gt_str_array_add(infiles, fas_path);
  gt_encseq_encoder_disable_description_support(encoder);
  gt_log_log("encode prefiltered reads to encseq %s", indexname);
  had_err = gt_encseq_encoder_encode(encoder, infiles, indexname, err);
  gt_encseq_encoder_delete(encoder);
  gt_str_delete(fas_path);
  gt_str_array_delete(infiles);
  return had_err;
}

static int gt_readjoiner_prefilter_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    void *tool_arguments, GtError *err)
{
  GtReadjoinerPrefilterArguments *arguments = tool_arguments;
  GtContfinder *contfinder;
  int had_err = 0;
  bool varlen;
  GtStr *cntlistfilename = NULL;
  GtStr *sepposfilename = NULL;
  unsigned long i;
  unsigned long input_nofreads, output_nofreads;
  GtLogger *default_logger, *verbose_logger;

  default_logger =
    gt_logger_new(!arguments->quiet, GT_LOGGER_DEFLT_PREFIX, stdout);
  verbose_logger =
    gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);

  gt_logger_log(default_logger, "gt readjoiner prefilter");

  if (gt_str_length(arguments->readset) == 0)
    gt_str_append_cstr(arguments->readset, gt_str_array_get(arguments->db, 0));
  gt_logger_log(verbose_logger, "readset name = %s",
      gt_str_get(arguments->readset));

  if (arguments->verbose)
  {
    GtStr *inputfileslist = gt_str_new();
    for (i = 0; i < gt_str_array_size(arguments->db); i++)
    {
      gt_str_append_cstr(inputfileslist, gt_str_array_get(arguments->db, i));
      if (i + 1 < gt_str_array_size(arguments->db))
        gt_str_append_cstr(inputfileslist, ", ");
    }
    gt_logger_log(verbose_logger, "input files = %s",
        gt_str_get(inputfileslist));
    gt_str_delete(inputfileslist);
  }

  if (arguments->cntlist)
  {
    cntlistfilename = gt_str_new_cstr(gt_str_get(arguments->readset));
    gt_str_append_cstr(cntlistfilename, GT_READJOINER_SUFFIX_CNTLIST);
  }
  if (arguments->seppos)
  {
    sepposfilename = gt_str_new_cstr(gt_str_get(arguments->readset));
    gt_str_append_cstr(sepposfilename, GT_READJOINER_SUFFIX_SEPPOS);
  }

  contfinder = gt_contfinder_new(arguments->db, arguments->readset,
      arguments->encseqall, err);
  varlen = gt_contfinder_read_length(contfinder) == 0;

  input_nofreads = gt_contfinder_nofseqs(contfinder) +
    gt_contfinder_nofdiscarded(contfinder);
  output_nofreads = input_nofreads;
  gt_logger_log(default_logger, "number of reads in complete readset = %lu",
      input_nofreads);
  if (varlen)
    gt_logger_log(verbose_logger, "read length = variable");
  else
    gt_logger_log(verbose_logger, "read length = %lu",
        gt_contfinder_read_length(contfinder));
  gt_logger_log(verbose_logger, "total length of complete readset = %lu",
      gt_contfinder_totallength_without_sep(contfinder) +
      gt_contfinder_discarded_length(contfinder));
  gt_logger_log(verbose_logger, "reads with ambiguities = %lu "
      "[%.2f %% of input]", gt_contfinder_nofdiscarded(contfinder),
      (float)gt_contfinder_nofdiscarded(contfinder) * 100 /
      (float)input_nofreads);
  if (!arguments->verbose)
    gt_logger_log(default_logger, "reads with ambiguities = %lu",
        gt_contfinder_nofdiscarded(contfinder));
  output_nofreads -= gt_contfinder_nofdiscarded(contfinder);

  if (contfinder == NULL)
    had_err = -1;
  else if (!arguments->encodeonly)
  {
    had_err = gt_contfinder_run(contfinder,
        !arguments->singlestrand, NULL,
        arguments->twobit ? GT_CONTFINDER_2BIT : (arguments->seqnums
          ? GT_CONTFINDER_SEQNUMS : (arguments->fasta ? GT_CONTFINDER_FASTA
            : GT_CONTFINDER_QUIET)), arguments->sorted, arguments->cntlist ?
        gt_str_get(cntlistfilename) : NULL, arguments->seppos ?
        gt_str_get(sepposfilename) : NULL, arguments->encseq, err);
    gt_logger_log(verbose_logger, "contained reads = %lu "
        "[%.2f %% of input]",
        gt_contfinder_nofcontained(contfinder),
        (float)gt_contfinder_nofcontained(contfinder) * 100 /
        (float)input_nofreads);
    if (!had_err)
    {
      if (!arguments->verbose)
      {
        gt_logger_log(default_logger, "contained reads = %lu",
            gt_contfinder_nofcontained(contfinder));
      }
      output_nofreads -= gt_contfinder_nofcontained(contfinder);
      gt_logger_log(default_logger, "number of reads in filtered readset = %lu",
                    output_nofreads);

      if (arguments->encseq)
      {
        gt_logger_log(verbose_logger, "encseq saved: %s.(esq|al1%s)",
                      gt_str_get(arguments->readset), varlen ? "|ssp" : "");
      }
      if (arguments->cntlist)
      {
        gt_logger_log(verbose_logger, "contained reads list saved: %s",
                      gt_str_get(cntlistfilename));
      }
      if (arguments->seppos)
      {
        gt_logger_log(verbose_logger, "separator positions saved: %s",
                      gt_str_get(sepposfilename));
      }
    }
  }
  if (!had_err)
  {
    if (arguments->testrs)
    {
      gt_contfinder_radixsort_eqlen_tester(contfinder);
    }
    gt_contfinder_delete(contfinder);
    if (arguments->encseq && varlen)
    {
      had_err
        = gt_prefilter_output_varlen_encseq(gt_str_get(arguments->readset),err);
    }
  }
  gt_str_delete(cntlistfilename);
  gt_str_delete(sepposfilename);
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
