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
#include "match/rdj-contigpaths.h"
#include "match/rdj-cntlist.h"
#include "match/rdj-spmlist.h"
#include "match/rdj-strgraph.h"
#include "match/rdj-filesuf-def.h"
#include "match/rdj-version.h"
#include "tools/gt_readjoiner_assembly.h"

typedef struct {
  bool verbose, quiet;
  unsigned int minmatchlength;
  unsigned int lengthcutoff, depthcutoff;
  GtStr  *readset, *buffersizearg;
  bool errors, paths2seq, redtrans, elendistri;
  unsigned int deadend, bubble, deadend_depth;
  GtOption *refoptionbuffersize;
  unsigned long buffersize;
  unsigned int nspmfiles;
} GtReadjoinerAssemblyArguments;

static void* gt_readjoiner_assembly_arguments_new(void)
{
  GtReadjoinerAssemblyArguments *arguments = gt_calloc((size_t)1,
      sizeof *arguments);
  arguments->readset = gt_str_new();
  arguments->buffersizearg = gt_str_new();
  arguments->buffersize = 0UL; /* in bytes */
  return arguments;
}

static void gt_readjoiner_assembly_arguments_delete(void *tool_arguments)
{
  GtReadjoinerAssemblyArguments *arguments = tool_arguments;
  if (!arguments)
    return;
  gt_str_delete(arguments->readset);
  gt_str_delete(arguments->buffersizearg);
  gt_option_delete(arguments->refoptionbuffersize);
  gt_free(arguments);
}

static GtOptionParser* gt_readjoiner_assembly_option_parser_new(
    void *tool_arguments)
{
  GtReadjoinerAssemblyArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *errors_option, *deadend_option, *v_option,
           *q_option, *bubble_option, *deadend_depth_option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...]",
      "Construct string graph and output contigs.");

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

  /* -depthcutoff */
  option = gt_option_new_uint_min("depthcutoff", "specify the minimal "
      "number of nodes in a contig",
      &arguments->depthcutoff, 3U, 1U);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -lengthcutoff */
  option = gt_option_new_uint_min("lengthcutoff", "specify the minimal "
      "length of a contig",
      &arguments->lengthcutoff, 100U, 1U);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -redtrans */
  option = gt_option_new_bool("redtrans", "reduce transitive edges",
      &arguments->redtrans, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -errors */
  errors_option = gt_option_new_bool("errors", "search graph features which "
      "may originate from sequencing errors and remove them",
      &arguments->errors, false);
  gt_option_is_extended_option(errors_option);
  gt_option_parser_add_option(op, errors_option);

  /* -bubble */
  bubble_option = gt_option_new_uint("bubble", "number of rounds of p-bubble "
      "removal to perform", &arguments->bubble, 3U);
  gt_option_is_extended_option(bubble_option);
  gt_option_imply(bubble_option, errors_option);
  gt_option_parser_add_option(op, bubble_option);

  /* -deadend */
  deadend_option = gt_option_new_uint("deadend", "number of rounds of "
      "dead end removal to perform a dead end",
      &arguments->deadend, 10U);
  gt_option_is_extended_option(deadend_option);
  gt_option_imply(deadend_option, errors_option);
  gt_option_parser_add_option(op, deadend_option);

  /* -deadend-depth */
  deadend_depth_option = gt_option_new_uint_min("deadend-depth", "specify the "
      "maximal depth of a path to an end-vertex by which the path shall be "
      "considered a dead end",
      &arguments->deadend_depth, 10U, 1U);
  gt_option_is_extended_option(deadend_depth_option);
  gt_option_imply(deadend_depth_option, errors_option);
  gt_option_parser_add_option(op, deadend_depth_option);

  /* -paths2seq */
  option = gt_option_new_bool("paths2seq", "read <indexname>"
      GT_READJOINER_SUFFIX_CONTIG_PATHS " and write "
      "<indexname>" GT_READJOINER_SUFFIX_CONTIGS,
      &arguments->paths2seq, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -buffersize */
  option = gt_option_new_string("buffersize", "specify size for read buffer"
      " of paths2seq phase (in bytes, the keywords 'MB' and 'GB' are allowed)",
                       arguments->buffersizearg, NULL);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);
  arguments->refoptionbuffersize = gt_option_ref(option);

  /* -v */
  v_option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, v_option);

  /* -q */
  q_option = gt_option_new_bool("q", "suppress standard output messages",
      &arguments->quiet, false);
  gt_option_parser_add_option(op, q_option);
  gt_option_exclude(q_option, v_option);

  /* -elendistri */
  option = gt_option_new_bool("elendistri", "output edges lenght"
      "distribution to <indexname>" GT_READJOINER_SUFFIX_ELEN_DISTRI,
      &arguments->elendistri, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_version_func(op, gt_readjoiner_show_version);
  gt_option_parser_set_max_args(op, 0);

  return op;
}

static int gt_readjoiner_assembly_arguments_check(GT_UNUSED int rest_argc,
    void *tool_arguments, GtError *err)
{
  GtReadjoinerAssemblyArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  if (had_err == 0 && gt_option_is_set(arguments->refoptionbuffersize))
    if (gt_option_parse_spacespec(&arguments->buffersize,
          "buffersize", arguments->buffersizearg, err) != 0)
      had_err = -1;
  return had_err;
}

#define GT_READJOINER_MSG_COUNTSPM \
  "calculate edges space for each vertex"
#define GT_READJOINER_MSG_BUILDSG \
  "build string graph"
#define GT_READJOINER_MSG_REDTRANS \
  "reduce transitive edges"
#define GT_READJOINER_MSG_CLEANSG \
  "correct sequencing errors"
#define GT_READJOINER_MSG_TRAVERSESG \
  "save contig paths"
#define GT_READJOINER_MSG_PUMPENCSEQ \
  "pump encseq through cache"
#define GT_READJOINER_MSG_OUTPUTCONTIGS \
  "save contig sequences"

static int gt_readjoiner_assembly_count_spm(const char *readset, bool eqlen,
    unsigned int minmatchlength, unsigned int nspmfiles, GtStrgraph *strgraph,
    GtBitsequence *contained, GtLogger *default_logger, GtError *err)
{
  GtSpmprocSkipData skipdata;
  int had_err = 0;
  unsigned int i;
  GtStr *filename = gt_str_new();
  gt_logger_log(default_logger, GT_READJOINER_MSG_COUNTSPM);
  if (!eqlen)
  {
    skipdata.out.e.proc = gt_spmproc_strgraph_count;
    skipdata.to_skip = contained;
    skipdata.out.e.data = strgraph;
  }
  for (i = 0; i < nspmfiles; i++)
  {
    gt_str_append_cstr(filename, readset);
    gt_str_append_char(filename, '.');
    gt_str_append_uint(filename, i);
    gt_str_append_cstr(filename, GT_READJOINER_SUFFIX_SPMLIST);
    had_err = gt_spmlist_parse(gt_str_get(filename),
        (unsigned long)minmatchlength,
        eqlen ? gt_spmproc_strgraph_count : gt_spmproc_skip,
        eqlen ? (void*)strgraph : (void*)&skipdata, err);
    gt_str_reset(filename);
  }
  gt_str_delete(filename);
  return had_err;
}

static int gt_readjoiner_assembly_error_correction(GtStrgraph *strgraph,
    unsigned int bubble, unsigned int deadend, unsigned int deadend_depth,
    GtLogger *verbose_logger)
{
  unsigned int i;
  unsigned long retval, retval_sum;
  gt_logger_log(verbose_logger, "remove p-bubbles");

  retval_sum = 0;
  retval = 1UL;
  for (i = 0; i < bubble && retval > 0; i++)
  {
    retval = gt_strgraph_redpbubbles(strgraph, 0, 1UL, false);
    retval_sum += retval;
    gt_logger_log(verbose_logger, "removed p-bubble edges [round %u] = %lu",
        i + 1, retval);
  }
  gt_logger_log(verbose_logger, "removed p-bubble edges [%u rounds] = %lu", i,
      retval_sum);
  gt_logger_log(verbose_logger, "remove dead-end paths");

  retval_sum = 0;
  retval = 1UL;
  for (i = 0; i < deadend && retval > 0; i++)
  {
    retval = gt_strgraph_reddepaths(strgraph, (unsigned long)deadend_depth,
        false);
    retval_sum += retval;
    gt_logger_log(verbose_logger, "removed dead-end path edges [round %u] = "
        "%lu", i + 1, retval);
  }
  gt_logger_log(verbose_logger,
      "removed dead-end path edges [%u rounds] = %lu", i, retval);
  return 0;
}

static void gt_readjoiner_assembly_pump_encseq_through_cache(
    const GtEncseq *encseq)
{
  const GtTwobitencoding *twobitencoding = gt_encseq_twobitencoding_export(
      encseq);
  uint64_t sum = 0; /* compute the sum, so that the compiler does no remove the
                       code accessing twobitencoding during optimization */
  unsigned long idx, totallength = gt_encseq_total_length(encseq),
                numofunits = ! gt_encseq_is_mirrored(encseq)
                  ? gt_unitsoftwobitencoding(totallength)
                  : gt_unitsoftwobitencoding((totallength - 1)/2);
  for (idx = 0; idx < numofunits; idx++)
    sum += twobitencoding[idx];
  gt_assert(sum > 0);
#ifndef S_SPLINT_S
  gt_log_log("encseq codes-sum: %"PRIu64, sum);
#endif
}

static int gt_readjoiner_assembly_paths2seq(const char *readset,
    unsigned long lengthcutoff, unsigned long buffersize,
    GtLogger *default_logger, GtTimer **timer, GtError *err)
{
  int had_err;
  GtEncseqLoader *el = gt_encseq_loader_new();
  GtEncseq *reads;

  if (gt_showtime_enabled())
  {
    gt_assert(timer != NULL);
    if (*timer == NULL) /* paths2seq */
    {
      *timer = gt_timer_new_with_progress_description(
          GT_READJOINER_MSG_PUMPENCSEQ);
      gt_timer_show_cpu_time_by_progress(*timer);
      gt_timer_start(*timer);
    }
    else
      gt_timer_show_progress(*timer, GT_READJOINER_MSG_PUMPENCSEQ, stdout);
  }
  gt_logger_log(default_logger, GT_READJOINER_MSG_PUMPENCSEQ);
  gt_encseq_loader_drop_description_support(el);
  gt_encseq_loader_disable_autosupport(el);
  gt_encseq_loader_mirror(el);
  reads = gt_encseq_loader_load(el, readset, err);
  gt_assert(reads != NULL);
  gt_readjoiner_assembly_pump_encseq_through_cache(reads);
  if (gt_showtime_enabled())
    gt_timer_show_progress(*timer, GT_READJOINER_MSG_OUTPUTCONTIGS, stdout);
  gt_logger_log(default_logger, GT_READJOINER_MSG_OUTPUTCONTIGS);
  had_err = gt_contigpaths_to_fasta(readset, GT_READJOINER_SUFFIX_CONTIG_PATHS,
      GT_READJOINER_SUFFIX_CONTIGS, reads, lengthcutoff, false,
      (size_t)buffersize, default_logger, err);
  gt_encseq_delete(reads);
  gt_encseq_loader_delete(el);
  return had_err;
}

static inline void gt_readjoiner_assembly_show_current_space(const char *label)
{
  unsigned long m, f;
  if (gt_ma_bookkeeping_enabled())
  {
    m = gt_ma_get_space_current();
    f = gt_fa_get_space_current();
    gt_log_log("used space %s: %.2f MB (ma: %.2f MB; fa: %.2f MB)",
        label == NULL ? "" : label, GT_MEGABYTES(m + f), GT_MEGABYTES(m),
        GT_MEGABYTES(f));
  }
}

static int gt_readjoiner_assembly_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    void *tool_arguments, GtError *err)
{
  GtReadjoinerAssemblyArguments *arguments = tool_arguments;
  GtLogger *verbose_logger, *default_logger;
  GtEncseqLoader *el;
  GtEncseq *reads;
  GtTimer *timer = NULL;
  GtStrgraph *strgraph = NULL;
  GtBitsequence *contained = NULL;
  const char *readset = gt_str_get(arguments->readset);
  bool eqlen;
  unsigned long nreads, tlen, rlen;
  int had_err = 0;

  gt_assert(arguments);
  gt_error_check(err);
  default_logger = gt_logger_new(!arguments->quiet, GT_LOGGER_DEFLT_PREFIX,
      stdout);
  gt_logger_log(default_logger,
      "gt readjoiner assembly (version "GT_READJOINER_VERSION")");
  verbose_logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX,
      stdout);
  gt_logger_log(verbose_logger, "verbose output activated");
  gt_logger_log(verbose_logger, "readset name = %s", readset);

  if (!arguments->paths2seq)
  {
    el = gt_encseq_loader_new();
    gt_encseq_loader_drop_description_support(el);
    gt_encseq_loader_disable_autosupport(el);
    reads = gt_encseq_loader_load(el, readset, err);
    gt_assert(reads != NULL);
    eqlen = gt_encseq_accesstype_get(reads) == GT_ACCESS_TYPE_EQUALLENGTH;
    nreads = gt_encseq_num_of_sequences(reads);
    gt_logger_log(default_logger, "number of reads in filtered readset = %lu",
        nreads);
    tlen = gt_encseq_total_length(reads) - nreads + 1;
    gt_logger_log(verbose_logger, "total length of filtered readset = %lu",
        tlen);

    if (eqlen)
    {
      rlen = gt_encseq_seqlength(reads, 0);
      gt_logger_log(verbose_logger, "read length = %lu", rlen);
      gt_encseq_delete(reads);
      reads = NULL;
    }
    else
    {
      unsigned int i;
      unsigned long nofreads;
      GtStr *filename = gt_str_clone(arguments->readset);
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
      rlen = 0;
      gt_logger_log(verbose_logger, "read length = variable");
      gt_assert(reads != NULL);
    }

    if (had_err == 0)
    {
      if (arguments->minmatchlength > 0)
        gt_logger_log(verbose_logger, "SPM length cutoff = %u",
            arguments->minmatchlength);

      strgraph = gt_strgraph_new(nreads);
      if (gt_showtime_enabled())
      {
        timer = gt_timer_new_with_progress_description(
            GT_READJOINER_MSG_COUNTSPM);
        gt_timer_start(timer);
        gt_timer_show_cpu_time_by_progress(timer);
      }
      had_err = gt_readjoiner_assembly_count_spm(readset, eqlen,
         arguments->minmatchlength, arguments->nspmfiles, strgraph, contained,
         default_logger, err);
      gt_readjoiner_assembly_show_current_space("(edges counted)");
      if (gt_showtime_enabled())
        gt_timer_show_progress(timer, GT_READJOINER_MSG_BUILDSG, stdout);
      gt_logger_log(default_logger, GT_READJOINER_MSG_BUILDSG);
    }
    if (had_err == 0)
    {
      gt_assert((eqlen && rlen > 0 && reads == NULL) ||
          (!eqlen && rlen == 0 && reads != NULL));
      gt_strgraph_allocate_graph(strgraph, rlen, reads);
      gt_readjoiner_assembly_show_current_space("(graph allocated)");
      had_err = gt_strgraph_load_spm_from_file(strgraph,
            (unsigned long)arguments->minmatchlength, arguments->redtrans,
            contained, readset, arguments->nspmfiles,
            GT_READJOINER_SUFFIX_SPMLIST, err);
    }
    if (had_err == 0)
    {
      if (arguments->elendistri)
        gt_strgraph_show_edge_lengths_distribution(strgraph, readset,
            GT_READJOINER_SUFFIX_ELEN_DISTRI);
      gt_strgraph_log_stats(strgraph, verbose_logger);
      gt_strgraph_log_space(strgraph);
    }

    if (!eqlen && reads != NULL && !arguments->errors)
    {
      gt_encseq_delete(reads);
      reads = NULL;
      if (had_err == 0)
        gt_strgraph_set_encseq(strgraph, NULL);
    }

    if (had_err == 0 && arguments->redtrans)
    {
      if (gt_showtime_enabled())
        gt_timer_show_progress(timer, GT_READJOINER_MSG_REDTRANS, stdout);
      gt_strgraph_sort_edges_by_len(strgraph, false);
      (void)gt_strgraph_redtrans(strgraph, false);
      (void)gt_strgraph_redself(strgraph, false);
      (void)gt_strgraph_redwithrc(strgraph, false);
      gt_strgraph_log_stats(strgraph, verbose_logger);
    }

    if (had_err == 0 && arguments->errors)
    {
      if (gt_showtime_enabled())
        gt_timer_show_progress(timer, GT_READJOINER_MSG_CLEANSG, stdout);
      gt_logger_log(default_logger, GT_READJOINER_MSG_CLEANSG);
      had_err = gt_readjoiner_assembly_error_correction(strgraph,
          arguments->bubble, arguments->deadend, arguments->deadend_depth,
          verbose_logger);
    }

    if (!eqlen && reads != NULL)
    {
      gt_encseq_delete(reads);
      reads = NULL;
      if (had_err == 0)
        gt_strgraph_set_encseq(strgraph, NULL);
    }

    if (had_err == 0)
    {
      if (gt_showtime_enabled())
        gt_timer_show_progress(timer, GT_READJOINER_MSG_TRAVERSESG, stdout);
      gt_logger_log(default_logger, GT_READJOINER_MSG_TRAVERSESG);
      gt_readjoiner_assembly_show_current_space("(before traversal)");
      gt_strgraph_spell(strgraph, (unsigned long)arguments->depthcutoff,
          (unsigned long)arguments->lengthcutoff, false, readset,
          GT_READJOINER_SUFFIX_CONTIG_PATHS, NULL, true, false,
          verbose_logger);
    }

    if (contained != NULL)
      gt_free(contained);
    gt_strgraph_delete(strgraph);
    strgraph = NULL;
    gt_assert(reads == NULL);
    gt_encseq_loader_delete(el);
  }

  if (had_err == 0)
  {
    gt_readjoiner_assembly_show_current_space("(before paths2seq)");
    had_err = gt_readjoiner_assembly_paths2seq(readset,
        (unsigned long)arguments->lengthcutoff,
        arguments->buffersize, default_logger, &timer, err);
  }

  if (gt_showtime_enabled())
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
  gt_logger_delete(default_logger);
  gt_logger_delete(verbose_logger);
  return had_err;
}

GtTool* gt_readjoiner_assembly(void)
{
  return gt_tool_new(gt_readjoiner_assembly_arguments_new,
                  gt_readjoiner_assembly_arguments_delete,
                  gt_readjoiner_assembly_option_parser_new,
                  gt_readjoiner_assembly_arguments_check,
                  gt_readjoiner_assembly_runner);
}
