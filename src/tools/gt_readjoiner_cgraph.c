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
#include "match/rdj-contigs-graph.h"
#include "match/rdj-filesuf-def.h"
#include "tools/gt_readjoiner_cgraph.h"

typedef struct {
  bool verbose, quiet, dot, smp1, smp2, dotdel, ext;
  GtStr  *readset;
  GtStrArray *subgraph;
  GtUword subgraph_depth;
} GtReadjoinerCgraphArguments;

static void* gt_readjoiner_cgraph_arguments_new(void)
{
  GtReadjoinerCgraphArguments *arguments = gt_calloc((size_t)1,
      sizeof *arguments);
  arguments->readset = gt_str_new();
  arguments->subgraph = gt_str_array_new();
  return arguments;
}

static void gt_readjoiner_cgraph_arguments_delete(void *tool_arguments)
{
  GtReadjoinerCgraphArguments *arguments = tool_arguments;
  if (!arguments)
    return;
  gt_str_delete(arguments->readset);
  gt_str_array_delete(arguments->subgraph);
  gt_free(arguments);
}

static GtOptionParser* gt_readjoiner_cgraph_option_parser_new(
    void *tool_arguments)
{
  GtReadjoinerCgraphArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *v_option, *subgraph_option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...]",
      "Construct contigs graph.");

  /* -readset */
  option = gt_option_new_string("readset", "specify the readset name",
      arguments->readset, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -v */
  v_option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, v_option);

  /* -q */
  option = gt_option_new_bool("q", "suppress standard output messages",
      &arguments->quiet, false);
  gt_option_exclude(option, v_option);
  gt_option_parser_add_option(op, option);

  /* -subgraph */
  subgraph_option = gt_option_new_string_array("subgraph",
      "output subgraph for the specified contigs "
      "in GraphViz format\n"
      "(<readset>"GT_READJOINER_SUFFIX_CG_SUB_DOT")",
      arguments->subgraph);
  gt_option_parser_add_option(op, subgraph_option);

  /* -dot */
  option = gt_option_new_bool("dot", "output complete graph in dot format",
      &arguments->dot, false);
  gt_option_parser_add_option(op, option);

  /* -smp1 */
  option = gt_option_new_bool("smp1", "simplify (algorithm 1)",
      &arguments->smp1, false);
  gt_option_parser_add_option(op, option);

  /* -smp2 */
  option = gt_option_new_bool("smp2", "simplify (algorithm 2)",
      &arguments->smp2, false);
  gt_option_parser_add_option(op, option);

  /* -ext */
  option = gt_option_new_bool("ext", "extend after simplify step",
      &arguments->ext, true);
  gt_option_parser_add_option(op, option);

  /* -dotdel */
  option = gt_option_new_bool("dotdel", "show deleted edges "
      "und vertices in dot output", &arguments->dotdel, false);
  gt_option_parser_add_option(op, option);

  /* -sdepth */
  option = gt_option_new_uword("sdepth", "subgraph depth\n"
      "(use 0 for infinite depth == whole connected component)",
      &(arguments->subgraph_depth), 1UL);
  gt_option_imply(option, subgraph_option);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_max_args(op, 0);

  return op;
}

#define GT_READJOINER_CGRAPH_GET_FP(FP, BASENAME, SUFFIX, MODE)\
  if (had_err == 0)\
  {\
    FP = gt_fa_fopen_with_suffix(BASENAME, SUFFIX, MODE, err);\
    if (FP == NULL)\
      had_err = -1;\
  }

static int gt_readjoiner_cgraph_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    void *tool_arguments, GtError *err)
{
  GtReadjoinerCgraphArguments *arguments = tool_arguments;
  GtLogger *verbose_logger, *default_logger;
  const char *readset = gt_str_get(arguments->readset);
  int had_err = 0;
  GtContigsGraph *cg = NULL;
  FILE *cjl_i_fp = NULL, *cjl_o_fp = NULL, *j_fp = NULL, *rlt_fp = NULL,
       *di_fp = NULL;

  gt_assert(arguments);
  gt_error_check(err);
  default_logger = gt_logger_new(!arguments->quiet, GT_LOGGER_DEFLT_PREFIX,
      stdout);
  gt_logger_log(default_logger, "gt readjoiner cgraph");
  verbose_logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX,
      stdout);
  gt_logger_log(verbose_logger, "verbose output activated");
  gt_logger_log(verbose_logger, "readset name = %s", readset);

  GT_READJOINER_CGRAPH_GET_FP(cjl_i_fp, readset,
      GT_READJOINER_SUFFIX_CJ_I_LINKS, "r");
  GT_READJOINER_CGRAPH_GET_FP(cjl_o_fp, readset,
      GT_READJOINER_SUFFIX_CJ_O_LINKS, "r");
  GT_READJOINER_CGRAPH_GET_FP(j_fp, readset,
      GT_READJOINER_SUFFIX_JUNCTIONS, "r");
  GT_READJOINER_CGRAPH_GET_FP(rlt_fp, readset,
      GT_READJOINER_SUFFIX_READSLIBRARYTABLE, "r");
  GT_READJOINER_CGRAPH_GET_FP(di_fp, readset,
      GT_READJOINER_SUFFIX_DEPTHINFO, "r");
  if (!had_err)
  {
    cg = gt_contigs_graph_new(cjl_i_fp, cjl_o_fp, j_fp, rlt_fp, di_fp, err);
    if (cg == NULL)
      had_err = -1;
  }
  if (cjl_i_fp != NULL)
    gt_fa_fclose(cjl_i_fp);
  if (cjl_o_fp != NULL)
    gt_fa_fclose(cjl_o_fp);
  if (j_fp != NULL)
    gt_fa_fclose(j_fp);
  if (rlt_fp != NULL)
    gt_fa_fclose(rlt_fp);
  if (di_fp != NULL)
    gt_fa_fclose(di_fp);
  if (!had_err && arguments->smp1)
  {
    gt_contigs_graph_simplify(cg, false);
    if (arguments->ext)
      gt_contigs_graph_extend_contigs(cg, true);
  }
  if (!had_err && arguments->smp2)
  {
    gt_contigs_graph_simplify(cg, true);
    if (arguments->ext)
      gt_contigs_graph_extend_contigs(cg, false);
  }
  if (!had_err && arguments->dotdel)
    gt_contigs_graph_enable_dot_show_deleted(cg);
  if (!had_err && arguments->dot)
  {
    FILE *dot_fp = NULL;
    GT_READJOINER_CGRAPH_GET_FP(dot_fp, readset,
        GT_READJOINER_SUFFIX_CG_DOT, "w");
    if (!had_err)
    {
      GtFile *dot_gtfile;
      dot_gtfile = gt_file_new_from_fileptr(dot_fp);
      gt_contigs_graph_show_dot(cg, dot_gtfile);
      gt_file_delete_without_handle(dot_gtfile);
    }
    if (dot_fp != NULL)
      gt_fa_fclose(dot_fp);
  }
  if (!had_err)
  {
    GtUword nofcnums = gt_str_array_size(arguments->subgraph);
    if (nofcnums > 0)
    {
      FILE *dot_fp = NULL;
      GtUword i, *cnums;
      const char *cnum_str;
      cnums = gt_malloc(sizeof (*cnums) * nofcnums);
      for (i = 0; i < nofcnums && !had_err; i++)
      {
        cnum_str = gt_str_array_get(arguments->subgraph, i);
        if (sscanf(cnum_str, ""GT_WU"", cnums + i) != 1)
        {
          gt_error_set(err, "argument of option -subgraph is invalid");
          had_err = -1;
        }
      }
      GT_READJOINER_CGRAPH_GET_FP(dot_fp, readset,
          GT_READJOINER_SUFFIX_CG_SUB_DOT, "w");
      if (!had_err)
      {
        GtFile *dot_gtfile;
        dot_gtfile = gt_file_new_from_fileptr(dot_fp);
        had_err = gt_contigs_graph_show_dot_subgraph(cg, dot_gtfile, cnums,
            nofcnums, arguments->subgraph_depth, err);
        gt_file_delete_without_handle(dot_gtfile);
      }
      if (dot_fp != NULL)
        gt_fa_fclose(dot_fp);
      gt_free(cnums);
    }
  }
  if (!had_err)
  {
    FILE *paths_fp = NULL;
    GT_READJOINER_CGRAPH_GET_FP(paths_fp, readset,
        GT_READJOINER_SUFFIX_CG_PATHS, "w");
    if (!had_err)
      gt_contigs_graph_output_paths(cg, paths_fp);
    if (paths_fp != NULL)
      gt_fa_fclose(paths_fp);
  }
  gt_contigs_graph_delete(cg);
  gt_logger_delete(default_logger);
  gt_logger_delete(verbose_logger);
  return had_err;
}

GtTool* gt_readjoiner_cgraph(void)
{
  return gt_tool_new(gt_readjoiner_cgraph_arguments_new,
                  gt_readjoiner_cgraph_arguments_delete,
                  gt_readjoiner_cgraph_option_parser_new,
                  NULL,
                  gt_readjoiner_cgraph_runner);
}
