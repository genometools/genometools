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
#include "tools/gt_readjoiner_graph.h"

typedef struct {
  GtStr  *readset;
  bool limits, verbose, quiet, errors, redtrans, eld,
       save, dot, mdot, adj, spm, subgraph_mono, subgraph_extend,
       firstonly, fromB, toB, fromE, toE;
  unsigned int deadend, bubble, deadend_depth;
  GtStrArray *subgraph, *subgraph_other;
  GtUword subgraph_depth, from, to, minlen, maxlen;
} GtReadjoinerGraphArguments;

static void* gt_readjoiner_graph_arguments_new(void)
{
  GtReadjoinerGraphArguments *arguments = gt_calloc((size_t)1,
      sizeof *arguments);
  arguments->readset = gt_str_new();
  arguments->subgraph = gt_str_array_new();
  arguments->subgraph_other = gt_str_array_new();
  return arguments;
}

static void gt_readjoiner_graph_arguments_delete(void *tool_arguments)
{
  GtReadjoinerGraphArguments *arguments = tool_arguments;
  if (!arguments)
    return;
  gt_str_delete(arguments->readset);
  gt_str_array_delete(arguments->subgraph);
  gt_str_array_delete(arguments->subgraph_other);
  gt_free(arguments);
}

static GtOptionParser* gt_readjoiner_graph_option_parser_new(
    void *tool_arguments)
{
  GtReadjoinerGraphArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *errors_option, *deadend_option, *v_option,
           *q_option, *bubble_option, *deadend_depth_option,
           *subgraph_option, *sother_option, *from_option,
           *to_option, *fromb_option, *tob_option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...]",
      "Collect information or further process the string graph.");

  /* -readset */
  option = gt_option_new_string("readset", "specify the readset name",
      arguments->readset, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -dot */
  option = gt_option_new_bool("dot",
      "output the string graph as bidirected graph "
      "in GraphViz format\n"
      "(<readset>"GT_READJOINER_SUFFIX_SG_BI_DOT")",
      &arguments->dot, false);
  gt_option_parser_add_option(op, option);

  /* -mdot */
  option = gt_option_new_bool("mdot",
      "output the string graph as monodirected graph "
      "(B/E vertices) in GraphViz format\n"
      "(<readset>"GT_READJOINER_SUFFIX_SG_MONO_DOT")",
      &arguments->mdot, false);
  gt_option_parser_add_option(op, option);

  /* -subgraph */
  subgraph_option = gt_option_new_string_array("subgraph",
      "output subgraph for the specified reads "
      "in GraphViz format\n"
      "(<readset>"GT_READJOINER_SUFFIX_SG_SUB_DOT")",
      arguments->subgraph);
  gt_option_parser_add_option(op, subgraph_option);

  /* -sdepth */
  option = gt_option_new_uword("sdepth", "subgraph depth\n"
      "(use 0 for infinite depth == whole connected component)",
      &(arguments->subgraph_depth), 1UL);
  gt_option_imply(option, subgraph_option);
  gt_option_parser_add_option(op, option);

  /* -sextend */
  option = gt_option_new_bool("sextend", "extend subgraph "
      "beyond depth along internal vertices",
      &arguments->subgraph_extend, false);
  gt_option_imply(option, subgraph_option);
  gt_option_parser_add_option(op, option);

  /* -smono */
  option = gt_option_new_bool("smono", "subgraph shall be "
      "output as monodirectional graph",
      &arguments->subgraph_mono, false);
  gt_option_imply(option, subgraph_option);
  gt_option_parser_add_option(op, option);

  /* -sother */
  sother_option = gt_option_new_string_array("sother",
      "also add specified reads to subgraph, but coloring "
      "them differently\n"
      "(<readset>"GT_READJOINER_SUFFIX_SG_SUB_DOT")",
      arguments->subgraph_other);
  gt_option_imply(sother_option, subgraph_option);
  gt_option_parser_add_option(op, sother_option);

  /* -spm */
  option = gt_option_new_bool("spm",
      "reconstruct SPM list from string graph\n"
      "(<readset>"GT_READJOINER_SUFFIX_SG_SPMLIST")",
      &arguments->spm, false);
  gt_option_parser_add_option(op, option);

  /* -adj */
  option = gt_option_new_bool("adj",
      "output graph adjacence list\n"
      "(<readset>"GT_READJOINER_SUFFIX_SG_ADJLIST")",
      &arguments->adj, false);
  gt_option_parser_add_option(op, option);

  /* -eld */
  option = gt_option_new_bool("eld", "output edges lenght"
      "distribution\n"
      "(<readset>"GT_READJOINER_SUFFIX_SG_ELEN_DISTRI")",
      &arguments->eld, false);
  gt_option_parser_add_option(op, option);

  /* -redtrans */
  option = gt_option_new_bool("redtrans", "reduce transitive edges",
      &arguments->redtrans, false);
  gt_option_parser_add_option(op, option);

  /* -errors */
  errors_option = gt_option_new_bool("errors",
      "run error correction algorithms",
      &arguments->errors, false);
  gt_option_parser_add_option(op, errors_option);

  /* -bubble */
  bubble_option = gt_option_new_uint("bubble", "number of rounds of p-bubble "
      "removal to perform", &arguments->bubble, 3U);
  gt_option_is_extended_option(bubble_option);
  gt_option_imply(bubble_option, errors_option);
  gt_option_parser_add_option(op, bubble_option);

  /* -from/-to */
  from_option = gt_option_new_uword("from", "find connecting path from->to\n"
      "requires -to; optional -minlen/-maxlen/-firstonly",
      &(arguments->from), 0);
  to_option = gt_option_new_uword("to", "find connecting path from->to\n"
      "requires -from; optional -minlen/-maxlen/-firstonly",
      &(arguments->to), 0);
  gt_option_imply(to_option, from_option);
  gt_option_imply(from_option, to_option);
  gt_option_parser_add_option(op, from_option);
  gt_option_parser_add_option(op, to_option);

  /* -fromB */
  fromb_option = gt_option_new_bool("fromB",
      "use only B vertex of origin read\n"
      "requires -from", &arguments->fromB, false);
  gt_option_imply(fromb_option, from_option);
  gt_option_parser_add_option(op, fromb_option);

  /* -fromE */
  option = gt_option_new_bool("fromE",
      "use only E vertex of origin read\n"
      "requires -from", &arguments->fromE, false);
  gt_option_imply(option, from_option);
  gt_option_exclude(option, fromb_option);
  gt_option_parser_add_option(op, option);

  /* -toB */
  tob_option = gt_option_new_bool("toB",
      "use only B vertex of origin read\n"
      "requires -to", &arguments->toB, false);
  gt_option_imply(tob_option, to_option);
  gt_option_parser_add_option(op, tob_option);

  /* -toE */
  option = gt_option_new_bool("toE",
      "use only E vertex of origin read\n"
      "requires -to", &arguments->toE, false);
  gt_option_imply(option, to_option);
  gt_option_exclude(option, tob_option);
  gt_option_parser_add_option(op, option);

  /* -minlen */
  option = gt_option_new_uword("minlen",
      "minimum string length of path from->to\n"
      "requires -from and -to",
      &(arguments->minlen), 0);
  gt_option_imply(option, from_option);
  gt_option_parser_add_option(op, option);

  /* -maxlen */
  option = gt_option_new_uword("maxlen",
      "maximum string length of path from->to\n"
      "requires -from and -to",
      &(arguments->maxlen), 1000);
  gt_option_imply(option, from_option);
  gt_option_parser_add_option(op, option);

  /* -firstonly */
  option = gt_option_new_bool("firstonly",
      "output only first path from->to",
      &arguments->firstonly, false);
  gt_option_parser_add_option(op, option);

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

  /* -save */
  option = gt_option_new_bool("save", "replace the string graph file "
      "after processing",
      &arguments->save, false);
  gt_option_parser_add_option(op, option);

  /* -limits */
  option = gt_option_new_bool("limits", "show compilation-specific limits",
      &arguments->limits, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -v */
  v_option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, v_option);

  /* -q */
  q_option = gt_option_new_bool("q", "suppress standard output messages",
      &arguments->quiet, false);
  gt_option_exclude(q_option, v_option);
  gt_option_parser_add_option(op, q_option);

  gt_option_parser_set_max_args(op, 0);

  return op;
}

#define GT_READJOINER_GRAPH_MSG_REDTRANS \
  "reduce transitive edges"
#define GT_READJOINER_GRAPH_MSG_CLEANSG \
  "correct sequencing errors"
#define GT_READJOINER_GRAPH_MSG_LOADSG \
  "load string graph from file"
#define GT_READJOINER_GRAPH_MSG_SAVESG \
  "save string graph to file"
#define GT_READJOINER_GRAPH_MSG_BDOTSG \
  "output bidirectional string graph in GraphViz format"
#define GT_READJOINER_GRAPH_MSG_MDOTSG \
  "output monodirectional string graph in GraphViz format"
#define GT_READJOINER_GRAPH_MSG_SDOTSG \
  "output specified subgraph in GraphViz format"
#define GT_READJOINER_GRAPH_MSG_ADJSG \
  "output graph adjacence list"
#define GT_READJOINER_GRAPH_MSG_SPMSG \
  "output SPM list from graph"

#define GT_READJOINER_GRAPH_ERR_SUBGRAPH_GENERIC(OPT) \
  "argument to "OPT" must be a list of read numbers"

#define GT_READJOINER_GRAPH_ERR_SUBGRAPH \
  GT_READJOINER_GRAPH_ERR_SUBGRAPH_GENERIC("subgraph")

#define GT_READJOINER_GRAPH_ERR_SOTHER \
  GT_READJOINER_GRAPH_ERR_SUBGRAPH_GENERIC("sother")

static int gt_readjoiner_graph_error_correction(GtStrgraph *strgraph,
    unsigned int bubble, unsigned int deadend, unsigned int deadend_depth,
    GtLogger *verbose_logger)
{
  unsigned int i;
  GtUword retval, retval_sum;
  gt_logger_log(verbose_logger, "remove p-bubbles");

  retval_sum = 0;
  retval = 1UL;
  for (i = 0; i < bubble && retval > 0; i++)
  {
    retval = gt_strgraph_redpbubbles(strgraph, 0, 1UL, false);
    retval_sum += retval;
    gt_logger_log(verbose_logger, "removed p-bubble edges [round %u] = "GT_WU"",
        i + 1, retval);
  }
  gt_logger_log(verbose_logger, "removed p-bubble edges [%u rounds] = "GT_WU"",
                i, retval_sum);
  gt_logger_log(verbose_logger, "remove dead-end paths");

  retval_sum = 0;
  retval = 1UL;
  for (i = 0; i < deadend && retval > 0; i++)
  {
    retval = gt_strgraph_reddepaths(strgraph, (GtUword)deadend_depth,
        false);
    retval_sum += retval;
    gt_logger_log(verbose_logger, "removed dead-end path edges [round %u] = "
        ""GT_WU"", i + 1, retval);
  }
  gt_logger_log(verbose_logger,
      "removed dead-end path edges [%u rounds] = "GT_WU"", i, retval_sum);
  return 0;
}

static inline void gt_readjoiner_graph_show_current_space(const char *label)
{
  GtUword m, f;
  if (gt_ma_bookkeeping_enabled())
  {
    m = gt_ma_get_space_current();
    f = gt_fa_get_space_current();
    gt_log_log("used space %s: %.2f MB (ma: %.2f MB; fa: %.2f MB)",
        label == NULL ? "" : label, GT_MEGABYTES(m + f), GT_MEGABYTES(m),
        GT_MEGABYTES(f));
  }
}

static void gt_readjoiner_graph_load_graph(GtStrgraph **strgraph,
    GtEncseq *reads, const char *readset, GtUword rlen,
    GtLogger *default_logger, GtTimer *timer)
{
  *strgraph = gt_strgraph_new_from_file(reads, rlen, readset,
      GT_READJOINER_SUFFIX_SG);

  if (gt_showtime_enabled())
    gt_timer_show_progress(timer, GT_READJOINER_GRAPH_MSG_LOADSG, stdout);
  gt_logger_log(default_logger, GT_READJOINER_GRAPH_MSG_LOADSG);

  gt_readjoiner_graph_show_current_space("(graph loaded)");
}

static int gt_readjoiner_graph_show_subgraph(GtStrgraph *strgraph,
    GtStrArray *subgraph, GtStrArray *subgraph_other,
    GtUword subgraph_depth, bool subgraph_mono, bool subgraph_extend,
    const char *readsetname, GtLogger *default_logger, GtTimer *timer,
    GtError *err)
{
  int had_err = 0;
  GtUword i, *rlist, rlistsize = gt_str_array_size(subgraph),
                *olist = NULL, olistsize = gt_str_array_size(subgraph_other);
  const char *num;
  gt_assert(rlistsize > 0);
  if (gt_showtime_enabled())
    gt_timer_show_progress(timer, GT_READJOINER_GRAPH_MSG_SDOTSG, stdout);
  gt_logger_log(default_logger, GT_READJOINER_GRAPH_MSG_SDOTSG);
  rlist = gt_malloc(sizeof (*rlist) * rlistsize);
  for (i = 0; i < rlistsize && !had_err; i++)
  {
    num = gt_str_array_get(subgraph, i);
    if (sscanf(num, ""GT_WU"", rlist + i) != 1)
    {
      gt_error_set(err, GT_READJOINER_GRAPH_ERR_SUBGRAPH);
      had_err = -1;
    }
  }
  if (olistsize > 0)
  {
    olist = gt_malloc(sizeof (*olist) * olistsize);
    for (i = 0; i < olistsize && !had_err; i++)
    {
      num = gt_str_array_get(subgraph_other, i);
      if (sscanf(num, ""GT_WU"", olist + i) != 1)
      {
        gt_error_set(err, GT_READJOINER_GRAPH_ERR_SOTHER);
        had_err = -1;
      }
    }
  }
  if (!had_err)
    had_err = gt_strgraph_show_context(strgraph, subgraph_mono ?
        GT_STRGRAPH_DOT : GT_STRGRAPH_DOT_BI, readsetname,
        GT_READJOINER_SUFFIX_SG_SUB_DOT, rlist, rlistsize,
        olist, olistsize, subgraph_depth == 0 ? ULONG_MAX : subgraph_depth,
        subgraph_extend, err);
  gt_free(rlist);
  gt_free(olist);
  return had_err;
}

static int gt_readjoiner_graph_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    void *tool_arguments, GtError *err)
{
  GtReadjoinerGraphArguments *arguments = tool_arguments;
  GtLogger *verbose_logger, *default_logger;
  GtEncseqLoader *el;
  GtEncseq *reads;
  GtTimer *timer = NULL;
  GtStrgraph *strgraph = NULL;
  const char *readset = gt_str_get(arguments->readset);
  bool eqlen;
  GtUword nreads, tlen, rlen;
  int had_err = 0;

  gt_assert(arguments);
  gt_error_check(err);

  default_logger = gt_logger_new(!arguments->quiet, GT_LOGGER_DEFLT_PREFIX,
      stdout);
  gt_logger_log(default_logger, "gt readjoiner graph");
  verbose_logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX,
      stdout);
  gt_logger_log(verbose_logger, "verbose output activated");
  gt_logger_log(verbose_logger, "readset name = %s", readset);

  if (arguments->limits)
    gt_strgraph_show_limits();

  if (gt_showtime_enabled())
  {
    timer = gt_timer_new_with_progress_description(
        GT_READJOINER_GRAPH_MSG_LOADSG);
    gt_timer_start(timer);
    gt_timer_show_cpu_time_by_progress(timer);
  }

  el = gt_encseq_loader_new();
  gt_encseq_loader_drop_description_support(el);
  gt_encseq_loader_disable_autosupport(el);
  reads = gt_encseq_loader_load(el, readset, err);
  gt_assert(reads != NULL);
  eqlen = gt_encseq_accesstype_get(reads) == GT_ACCESS_TYPE_EQUALLENGTH;
  nreads = gt_encseq_num_of_sequences(reads);
  gt_logger_log(default_logger, "number of reads in filtered readset = "GT_WU"",
      nreads);
  tlen = gt_encseq_total_length(reads) - nreads + 1;
  gt_logger_log(verbose_logger, "total length of filtered readset = "GT_WU"",
      tlen);

  if (eqlen)
  {
    rlen = gt_encseq_seqlength(reads, 0);
    gt_logger_log(verbose_logger, "read length = "GT_WU"", rlen);
    gt_encseq_delete(reads);
    reads = NULL;
  }
  else
  {
    rlen = 0;
    gt_logger_log(verbose_logger, "read length = variable");
    gt_assert(reads != NULL);
  }

  if (had_err == 0)
  {
    gt_readjoiner_graph_load_graph(&strgraph, reads, readset, rlen,
        default_logger, timer);
  }

  if (had_err == 0)
  {
    if (arguments->eld)
      gt_strgraph_show_edge_lengths_distribution(strgraph, readset,
          GT_READJOINER_SUFFIX_SG_ELEN_DISTRI);
    gt_strgraph_log_stats(strgraph, verbose_logger);
    gt_strgraph_log_space(strgraph);
  }

  if (had_err == 0 && arguments->from)
  {
    GtStrgraphVtype from_vt, to_vt;
    gt_assert(eqlen || (reads != NULL));
    gt_encseq_mirror(reads, err);
    if (arguments->fromB)
      from_vt = GT_STRGRAPH_VTYPE_B;
    else if (arguments->fromE)
      from_vt = GT_STRGRAPH_VTYPE_E;
    else
      from_vt = GT_STRGRAPH_VTYPE_A;
    if (arguments->toB)
      to_vt = GT_STRGRAPH_VTYPE_B;
    else if (arguments->toE)
      to_vt = GT_STRGRAPH_VTYPE_E;
    else
      to_vt = GT_STRGRAPH_VTYPE_A;
    had_err = gt_strgraph_find_connecting_path(strgraph, arguments->from,
        from_vt, arguments->to, to_vt, arguments->minlen, arguments->maxlen,
        arguments->firstonly, readset, GT_READJOINER_SUFFIX_CONNECTING_PATHS,
        verbose_logger, err);
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
      gt_timer_show_progress(timer, GT_READJOINER_GRAPH_MSG_REDTRANS, stdout);
    gt_strgraph_sort_edges_by_len(strgraph, false);
    (void)gt_strgraph_redtrans(strgraph, false);
    (void)gt_strgraph_redself(strgraph, false);
    (void)gt_strgraph_redwithrc(strgraph, false);
    gt_strgraph_log_stats(strgraph, verbose_logger);
  }

  if (had_err == 0 && arguments->errors)
  {
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, GT_READJOINER_GRAPH_MSG_CLEANSG, stdout);
    gt_logger_log(default_logger, GT_READJOINER_GRAPH_MSG_CLEANSG);
    had_err = gt_readjoiner_graph_error_correction(strgraph,
        arguments->bubble, arguments->deadend, arguments->deadend_depth,
        verbose_logger);
  }

  if (had_err == 0 && arguments->save)
  {
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, GT_READJOINER_GRAPH_MSG_SAVESG, stdout);
    gt_logger_log(default_logger, GT_READJOINER_GRAPH_MSG_SAVESG);
    gt_strgraph_show(strgraph, GT_STRGRAPH_BIN,
        gt_str_get(arguments->readset), GT_READJOINER_SUFFIX_SG, false);
  }

  if (had_err == 0 && arguments->dot)
  {
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, GT_READJOINER_GRAPH_MSG_BDOTSG, stdout);
    gt_logger_log(default_logger, GT_READJOINER_GRAPH_MSG_BDOTSG);
    gt_strgraph_show(strgraph, GT_STRGRAPH_DOT_BI,
        gt_str_get(arguments->readset), GT_READJOINER_SUFFIX_SG_BI_DOT, false);
  }

  if (had_err == 0 && arguments->mdot)
  {
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, GT_READJOINER_GRAPH_MSG_MDOTSG, stdout);
    gt_logger_log(default_logger, GT_READJOINER_GRAPH_MSG_MDOTSG);
    gt_strgraph_show(strgraph, GT_STRGRAPH_DOT,
        gt_str_get(arguments->readset), GT_READJOINER_SUFFIX_SG_MONO_DOT,
        false);
  }

  if (had_err == 0 && gt_str_array_size(arguments->subgraph) > 0)
  {
    had_err = gt_readjoiner_graph_show_subgraph(strgraph, arguments->subgraph,
        arguments->subgraph_other, arguments->subgraph_depth,
        arguments->subgraph_mono, arguments->subgraph_extend,
        gt_str_get(arguments->readset), default_logger, timer, err);
  }

  if (had_err == 0 && arguments->adj)
  {
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, GT_READJOINER_GRAPH_MSG_ADJSG, stdout);
    gt_logger_log(default_logger, GT_READJOINER_GRAPH_MSG_ADJSG);
    gt_strgraph_show(strgraph, GT_STRGRAPH_ADJLIST,
        gt_str_get(arguments->readset), GT_READJOINER_SUFFIX_SG_ADJLIST, false);
  }

  if (had_err == 0 && arguments->spm)
  {
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, GT_READJOINER_GRAPH_MSG_SPMSG, stdout);
    gt_logger_log(default_logger, GT_READJOINER_GRAPH_MSG_SPMSG);
    gt_strgraph_show(strgraph, GT_STRGRAPH_SPM,
        gt_str_get(arguments->readset), GT_READJOINER_SUFFIX_SG_SPMLIST, false);
  }

  if (!eqlen && reads != NULL)
  {
    gt_encseq_delete(reads);
    reads = NULL;
    if (had_err == 0)
      gt_strgraph_set_encseq(strgraph, NULL);
  }

  gt_strgraph_delete(strgraph);
  strgraph = NULL;
  gt_assert(reads == NULL);
  gt_encseq_loader_delete(el);

  if (had_err == 0)
  {
    gt_readjoiner_graph_show_current_space("(final)");
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

GtTool* gt_readjoiner_graph(void)
{
  return gt_tool_new(gt_readjoiner_graph_arguments_new,
                  gt_readjoiner_graph_arguments_delete,
                  gt_readjoiner_graph_option_parser_new,
                  NULL,
                  gt_readjoiner_graph_runner);
}
