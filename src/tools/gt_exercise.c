/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/cstr_array.h"
#include "core/option.h"
#include "core/versionfunc.h"
#include "extended/toolbox.h"
#include "tools/gt_affinealign.h"
#include "tools/gt_align.h"
#include "tools/gt_assemblegreedy.h"
#include "tools/gt_blastenv.h"
#include "tools/gt_casino.h"
#include "tools/gt_coin.h"
#include "tools/gt_consensus_sa.h"
#include "tools/gt_exercise.h"
#include "tools/gt_fastaparser.h"
#include "tools/gt_linearalign.h"
#include "tools/gt_markovchain.h"
#include "tools/gt_matchcount.h"
#include "tools/gt_msaparse.h"
#include "tools/gt_multilcp.h"
#include "tools/gt_multiset_matching.h"
#include "tools/gt_neighborjoining.h"
#include "tools/gt_nussinov_rna_fold.h"
#include "tools/gt_qgramdist.h"
#include "tools/gt_scorefasta.h"
#include "tools/gt_scorematrix.h"
#include "tools/gt_swalign.h"
#include "tools/gt_translate.h"
#include "tools/gt_upgma.h"

static void* gt_exercise_arguments_new(void)
{
  GtToolbox *exercise_toolbox = gt_toolbox_new();
  gt_toolbox_add_tool(exercise_toolbox, "affinealign", gt_affinealign_tool());
  gt_toolbox_add_tool(exercise_toolbox, "align", gt_align());
  gt_toolbox_add_tool(exercise_toolbox, "assemblegreedy", gt_assemblegreedy());
  gt_toolbox_add_tool(exercise_toolbox, "blastenv", gt_blastenv());
  gt_toolbox_add(exercise_toolbox, "casino", gt_casino);
  gt_toolbox_add(exercise_toolbox, "coin", gt_coin);
  gt_toolbox_add_tool(exercise_toolbox, "consensus_sa", gt_consensus_sa_tool());
  gt_toolbox_add_tool(exercise_toolbox, "fastaparser", gt_fastaparser());
  gt_toolbox_add_tool(exercise_toolbox, "markovchain", gt_markovchain());
  gt_toolbox_add(exercise_toolbox, "matchcount", gt_matchcount);
  gt_toolbox_add(exercise_toolbox, "msaparse", gt_msaparse);
  gt_toolbox_add(exercise_toolbox, "msmatch", gt_multiset_matching);
  gt_toolbox_add(exercise_toolbox, "multilcp", gt_multilcp);
  gt_toolbox_add_tool(exercise_toolbox, "linearalign", gt_linearalign());
  gt_toolbox_add(exercise_toolbox, "neighborjoining", gt_neighborjoining);
  gt_toolbox_add(exercise_toolbox, "nussinov_rna_fold", gt_nussinov_rna_fold);
  gt_toolbox_add(exercise_toolbox, "qgramdist", gt_qgramdist);
  gt_toolbox_add_tool(exercise_toolbox, "scorefasta", gt_scorefasta());
  gt_toolbox_add(exercise_toolbox, "scorematrix", gt_scorematrix);
  gt_toolbox_add_tool(exercise_toolbox, "swalign", gt_swalign());
  gt_toolbox_add(exercise_toolbox, "translate", gt_translate);
  gt_toolbox_add(exercise_toolbox, "upgma", gt_upgma);
  return exercise_toolbox;
}

static void gt_exercise_arguments_delete(void *tool_arguments)
{
  GtToolbox *exercise_toolbox = tool_arguments;
  if (!exercise_toolbox) return;
  gt_toolbox_delete(exercise_toolbox);
}

static GtOptionParser* gt_exercise_option_parser_new(void *tool_arguments)
{
  GtToolbox *exercise_toolbox = tool_arguments;
  GtOptionParser *op;
  assert(exercise_toolbox);
  op = gt_option_parser_new("[option ...] exercise_tool_name [argument ...]",
                         "Call exercise tool with name exercise_tool_name and "
                         "pass argument(s) to it.");
  gt_option_parser_set_comment_func(op, gt_toolbox_show, exercise_toolbox);
  gt_option_parser_set_min_args(op, 1);
  return op;
}

static int gt_exercise_runner(int argc, const char **argv, int parsed_args,
                             void *tool_arguments, GtError *err)
{
  GtToolbox *exercise_toolbox = tool_arguments;
  GtToolfunc toolfunc;
  GtTool *tool = NULL;
  int had_err = 0;
  char **nargv = NULL;

  gt_error_check(err);
  assert(exercise_toolbox);

  /* get exercise */
  if (!gt_toolbox_has_tool(exercise_toolbox, argv[parsed_args])) {
    gt_error_set(err, "exercise '%s' not found; option -help lists possible "
                   "tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call exercise */
  if (!had_err) {
    if (!(toolfunc = gt_toolbox_get(exercise_toolbox, argv[parsed_args]))) {
      tool = gt_toolbox_get_tool(exercise_toolbox, argv[parsed_args]);
      assert(tool);
    }
    nargv = gt_cstr_array_prefix_first(argv + parsed_args,
                                       gt_error_get_progname(err));
    gt_error_set_progname(err, nargv[0]);
    if (toolfunc)
      had_err = toolfunc(argc - parsed_args, (const char**) nargv, err);
    else
      had_err = gt_tool_run(tool, argc - parsed_args, (const char**) nargv,
                            err);
  }

  /* free */
  gt_cstr_array_delete(nargv);

  return had_err;
}

GtTool* gt_exercise(void)
{
  return gt_tool_new(gt_exercise_arguments_new,
                  gt_exercise_arguments_delete,
                  gt_exercise_option_parser_new,
                  NULL,
                  gt_exercise_runner);
}
