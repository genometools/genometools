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

#include "core/encseq.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/fa.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/showtime.h"
#include "core/spacecalc.h"
#include "match/asqg_writer.h"
#include "match/rdj-contigpaths.h"
#include "match/rdj-cntlist.h"
#include "match/rdj-spmlist.h"
#include "match/rdj-strgraph.h"
#include "match/rdj-filesuf-def.h"
#include "match/rdj-version.h"
#include "tools/gt_readjoiner_asqg.h"

typedef struct {
  bool verbose, quiet;
  unsigned int minmatchlength;
  GtStr  *readset;
  bool gz, sg;
  unsigned int nspmfiles;
} GtReadjoinerAsqgArguments;

static void* gt_readjoiner_asqg_arguments_new(void)
{
  GtReadjoinerAsqgArguments *arguments = gt_calloc((size_t)1,
      sizeof *arguments);
  arguments->readset = gt_str_new();
  return arguments;
}

static void gt_readjoiner_asqg_arguments_delete(void *tool_arguments)
{
  GtReadjoinerAsqgArguments *arguments = tool_arguments;
  if (!arguments)
    return;
  gt_str_delete(arguments->readset);
  gt_free(arguments);
}

static GtOptionParser* gt_readjoiner_asqg_option_parser_new(
    void *tool_arguments)
{
  GtReadjoinerAsqgArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *v_option;
  gt_assert(arguments != NULL);

  /* init */
  op = gt_option_parser_new("[option ...]",
      "Output string graph in SGA asqg format.");

  /* -readset */
  option = gt_option_new_string("readset", "specify the readset name",
      arguments->readset, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -spmfiles */
  option = gt_option_new_uint_min("spmfiles", "number of SPM files to read\n"
      "this must be equal to the value of -j for the overlap phase",
      &arguments->nspmfiles, 1U, 1U);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -l */
  option = gt_option_new_uint_min("l", "specify the minimum SPM length",
      &arguments->minmatchlength, 0, 2U);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -gz */
  option = gt_option_new_bool("gz", "output gzipped file",
      &arguments->gz, false);
  gt_option_parser_add_option(op, option);

  /* -sg */
  option = gt_option_new_bool("sg", "first construct a Readjoiner string "
      "graph, then convert it into SGA format",
      &arguments->sg, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -v */
  v_option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, v_option);

  /* -q */
  option = gt_option_new_bool("q", "suppress standard output messages",
      &arguments->quiet, false);
  gt_option_exclude(option, v_option);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_version_func(op, gt_readjoiner_show_version);
  gt_option_parser_set_max_args(op, 0);

  return op;
}

#define GT_READJOINER_ASQG_MSG_INIT \
  "initialization"
#define GT_READJOINER_ASQG_MSG_CNT \
  "parse lists of contained reads"

/* encseq + spm -> asqg */
#define GT_READJOINER_ASQG_MSG_VERTICES \
  "output vertices in asqg format"
#define GT_READJOINER_ASQG_MSG_EDGES \
  "output edges in asqg format"

/* encseq + sg -> asqg */
#define GT_READJOINER_ASQG_MSG_COUNT \
  "build string graph (counting phase)"
#define GT_READJOINER_ASQG_MSG_INSERT \
  "build string graph (insertion phase)"
#define GT_READJOINER_ASQG_MSG_OUTPUT \
  "output string graph in asqg format"

static int gt_readjoiner_asqg_use_spmfiles(GtSpmproc proc, void *procdata,
    const char *readset, unsigned int minmatchlength, unsigned int nspmfiles,
    GtBitsequence *contained, GtError *err)
{
  int had_err = 0;
  GtSpmprocSkipData skipdata;
  unsigned int i;
  GtStr *filename = gt_str_new();
  if (contained != NULL)
  {
    skipdata.to_skip = contained;
    skipdata.out.e.proc = proc;
    skipdata.out.e.data = procdata;
  }
  for (i = 0; i < nspmfiles; i++)
  {
    gt_str_append_cstr(filename, readset);
    gt_str_append_char(filename, '.');
    gt_str_append_uint(filename, i);
    gt_str_append_cstr(filename, GT_READJOINER_SUFFIX_SPMLIST);
    had_err = gt_spmlist_parse(gt_str_get(filename),
        (unsigned long)minmatchlength,
        contained == NULL ? proc : gt_spmproc_skip,
        contained == NULL ? (void*)procdata : (void*)&skipdata, err);
    gt_str_reset(filename);
  }
  gt_str_delete(filename);
  return had_err;
}

static inline void gt_readjoiner_asqg_show_current_space(const char *label)
{
  unsigned long m, f;
  if (gt_ma_bookkeeping_enabled())
  {
    m = gt_ma_get_space_current();
    f = gt_fa_get_space_current();
    gt_log_log("used space after %s: %.2f MB (ma: %.2f MB; fa: %.2f MB)",
        label == NULL ? "" : label, GT_MEGABYTES(m + f), GT_MEGABYTES(m),
        GT_MEGABYTES(f));
  }
}

static int gt_readjoiner_asqg_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    void *tool_arguments, GtError *err)
{
  GtReadjoinerAsqgArguments *arguments = tool_arguments;
  GtLogger *verbose_logger, *default_logger;
  GtEncseqLoader *el;
  GtEncseq *reads;
  GtTimer *timer = NULL;
  GtStrgraph *strgraph = NULL;
  GtBitsequence *contained = NULL;
  const char *readset = gt_str_get(arguments->readset);
  bool eqlen;
  unsigned long nreads;
  int had_err = 0;

  if (gt_showtime_enabled())
  {
    timer = gt_timer_new_with_progress_description(
        GT_READJOINER_ASQG_MSG_INIT);
    gt_timer_start(timer);
    gt_timer_show_cpu_time_by_progress(timer);
  }
  gt_assert(arguments);
  gt_error_check(err);
  default_logger = gt_logger_new(!arguments->quiet, GT_LOGGER_DEFLT_PREFIX,
      stdout);
  gt_logger_log(default_logger, "gt readjoiner asqg");
  verbose_logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX,
      stdout);
  gt_logger_log(verbose_logger, "verbose output activated");
  gt_logger_log(verbose_logger, "readset name = %s", readset);
  if (arguments->minmatchlength > 0)
    gt_logger_log(verbose_logger, "SPM length cutoff = %u",
        arguments->minmatchlength);

  el = gt_encseq_loader_new();
  gt_encseq_loader_drop_description_support(el);
  gt_encseq_loader_disable_autosupport(el);
  reads = gt_encseq_loader_load(el, readset, err);
  nreads = gt_encseq_num_of_sequences(reads);
  gt_assert(reads != NULL);
  eqlen = gt_encseq_accesstype_get(reads) == GT_ACCESS_TYPE_EQUALLENGTH;
  gt_readjoiner_asqg_show_current_space(GT_READJOINER_ASQG_MSG_INIT);
  if (!eqlen)
  {
    unsigned int i;
    unsigned long nofreads;
    GtStr *filename = gt_str_clone(arguments->readset);
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, GT_READJOINER_ASQG_MSG_CNT, stdout);
    gt_logger_log(default_logger, GT_READJOINER_ASQG_MSG_CNT);
    gt_str_append_cstr(filename, ".0" GT_READJOINER_SUFFIX_CNTLIST);
    had_err = gt_cntlist_parse(gt_str_get(filename), true, &contained,
        &nofreads, err);
    for (i = 1U; i < arguments->nspmfiles && had_err == 0; i++)
    {
      unsigned long nofreads_i;
      gt_str_reset(filename);
      gt_str_append_str(filename, arguments->readset);
      gt_str_append_char(filename, '.');
      gt_str_append_uint(filename, i);
      gt_str_append_cstr(filename, GT_READJOINER_SUFFIX_CNTLIST);
      had_err = gt_cntlist_parse(gt_str_get(filename), false, &contained,
          &nofreads_i, err);
      gt_assert(nofreads == nofreads_i);
    }
    gt_str_delete(filename);
    gt_readjoiner_asqg_show_current_space(GT_READJOINER_ASQG_MSG_CNT);
  }
  if (!arguments->sg)
  {
    GtStr *filename = NULL;
    GtFile *file = NULL;
    GtAsqgWriter *aw = NULL;
    if (had_err == 0)
    {
      filename = gt_str_clone(arguments->readset);
      gt_str_append_cstr(filename, arguments->gz ? ".asqg.gz" : ".asqg");
      file = gt_file_open(arguments->gz ? GT_FILE_MODE_GZIP :
          GT_FILE_MODE_UNCOMPRESSED, gt_str_get(filename), "w", err);
      if (file == NULL)
        had_err = -1;
    }
    if (had_err == 0)
    {
      aw = gt_asqg_writer_new(file, reads);
      if (gt_showtime_enabled())
        gt_timer_show_progress(timer, GT_READJOINER_ASQG_MSG_VERTICES,
            stdout);
      gt_logger_log(default_logger, GT_READJOINER_ASQG_MSG_VERTICES);
      had_err = gt_asqg_writer_show_header(aw, 0.0,
          (unsigned long)arguments->minmatchlength,
          gt_str_get(arguments->readset), false, false, err);
    }
    if (had_err == 0)
    {
      had_err = gt_asqg_writer_show_vertices(aw, err);
      gt_readjoiner_asqg_show_current_space(GT_READJOINER_ASQG_MSG_VERTICES);
    }
    if (had_err == 0)
    {
      if (gt_showtime_enabled())
        gt_timer_show_progress(timer, GT_READJOINER_ASQG_MSG_EDGES, stdout);
      gt_logger_log(default_logger, GT_READJOINER_ASQG_MSG_EDGES);
      had_err = gt_readjoiner_asqg_use_spmfiles(gt_spmproc_show_asgq,
          aw, readset, arguments->minmatchlength, arguments->nspmfiles,
          contained, err);
      gt_readjoiner_asqg_show_current_space(GT_READJOINER_ASQG_MSG_EDGES);
    }
    gt_str_delete(filename);
    gt_file_delete(file);
    gt_asqg_writer_delete(aw);
  }
  else
  {
    if (had_err == 0)
    {
      if (gt_showtime_enabled())
        gt_timer_show_progress(timer, GT_READJOINER_ASQG_MSG_COUNT, stdout);
      gt_logger_log(default_logger, GT_READJOINER_ASQG_MSG_COUNT);
      strgraph = gt_strgraph_new(nreads);
      had_err = gt_readjoiner_asqg_use_spmfiles(gt_spmproc_strgraph_count,
          strgraph, readset, arguments->minmatchlength, arguments->nspmfiles,
          contained, err);
      gt_readjoiner_asqg_show_current_space(GT_READJOINER_ASQG_MSG_COUNT);
    }
    if (had_err == 0)
    {
      if (gt_showtime_enabled())
        gt_timer_show_progress(timer, GT_READJOINER_ASQG_MSG_INSERT, stdout);
      gt_logger_log(default_logger, GT_READJOINER_ASQG_MSG_INSERT);
      gt_strgraph_allocate_graph(strgraph,
          eqlen ? gt_encseq_seqlength(reads, 0) : 0,
          eqlen ? NULL : reads);
      had_err = gt_strgraph_load_spm_from_file(strgraph,
          (unsigned long)arguments->minmatchlength, false,
          contained, readset, arguments->nspmfiles,
          GT_READJOINER_SUFFIX_SPMLIST, err);
      gt_readjoiner_asqg_show_current_space(GT_READJOINER_ASQG_MSG_INSERT);
    }
    if (had_err == 0)
    {
      if (gt_showtime_enabled())
        gt_timer_show_progress(timer, GT_READJOINER_ASQG_MSG_OUTPUT, stdout);
      gt_logger_log(default_logger, GT_READJOINER_ASQG_MSG_OUTPUT);
      gt_strgraph_set_encseq(strgraph, reads);
      gt_strgraph_show(strgraph, arguments->gz ? GT_STRGRAPH_ASQG_GZ :
          GT_STRGRAPH_ASQG, gt_str_get(arguments->readset),
          arguments->gz ? ".asqg.gz" : ".asqg", false);
      gt_readjoiner_asqg_show_current_space(GT_READJOINER_ASQG_MSG_OUTPUT);
    }
  }
  if (gt_showtime_enabled())
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
  gt_free(contained);
  gt_strgraph_delete(strgraph);
  gt_encseq_loader_delete(el);
  gt_encseq_delete(reads);
  gt_logger_delete(default_logger);
  gt_logger_delete(verbose_logger);
  return had_err;
}

GtTool* gt_readjoiner_asqg(void)
{
  return gt_tool_new(gt_readjoiner_asqg_arguments_new,
                  gt_readjoiner_asqg_arguments_delete,
                  gt_readjoiner_asqg_option_parser_new,
                  NULL, gt_readjoiner_asqg_runner);
}
