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

#include "libgtcore/cstr_array.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/toolbox.h"
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
  Toolbox *exercise_toolbox = toolbox_new();
  toolbox_add_tool(exercise_toolbox, "affinealign", gt_affinealign());
  toolbox_add_tool(exercise_toolbox, "align", gt_align());
  toolbox_add_tool(exercise_toolbox, "assemblegreedy", gt_assemblegreedy());
  toolbox_add_tool(exercise_toolbox, "blastenv", gt_blastenv());
  toolbox_add(exercise_toolbox, "casino", gt_casino);
  toolbox_add(exercise_toolbox, "coin", gt_coin);
  toolbox_add_tool(exercise_toolbox, "consensus_sa", gt_consensus_sa());
  toolbox_add_tool(exercise_toolbox, "fastaparser", gt_fastaparser());
  toolbox_add_tool(exercise_toolbox, "markovchain", gt_markovchain());
  toolbox_add(exercise_toolbox, "matchcount", gt_matchcount);
  toolbox_add(exercise_toolbox, "msaparse", gt_msaparse);
  toolbox_add(exercise_toolbox, "msmatch", gt_multiset_matching);
  toolbox_add(exercise_toolbox, "multilcp", gt_multilcp);
  toolbox_add_tool(exercise_toolbox, "linearalign", gt_linearalign());
  toolbox_add(exercise_toolbox, "neighborjoining", gt_neighborjoining);
  toolbox_add(exercise_toolbox, "nussinov_rna_fold", gt_nussinov_rna_fold);
  toolbox_add(exercise_toolbox, "qgramdist", gt_qgramdist);
  toolbox_add_tool(exercise_toolbox, "scorefasta", gt_scorefasta());
  toolbox_add(exercise_toolbox, "scorematrix", gt_scorematrix);
  toolbox_add_tool(exercise_toolbox, "swalign", gt_swalign());
  toolbox_add(exercise_toolbox, "translate", gt_translate);
  toolbox_add(exercise_toolbox, "upgma", gt_upgma);
  return exercise_toolbox;
}

static void gt_exercise_arguments_delete(void *tool_arguments)
{
  Toolbox *exercise_toolbox = tool_arguments;
  if (!exercise_toolbox) return;
  toolbox_delete(exercise_toolbox);
}

static OptionParser* gt_exercise_option_parser_new(void *tool_arguments)
{
  Toolbox *exercise_toolbox = tool_arguments;
  OptionParser *op;
  assert(exercise_toolbox);
  op = option_parser_new("[option ...] exercise_tool_name [argument ...]",
                         "Call exercise tool with name exercise_tool_name and "
                         "pass argument(s) to it.");
  option_parser_set_comment_func(op, toolbox_show, exercise_toolbox);
  option_parser_set_min_args(op, 1);
  return op;
}

static int gt_exercise_runner(int argc, const char **argv, int parsed_args,
                             void *tool_arguments, Error *err)
{
  Toolbox *exercise_toolbox = tool_arguments;
  Toolfunc toolfunc;
  Tool *tool = NULL;
  int had_err = 0;
  char **nargv = NULL;

  error_check(err);
  assert(exercise_toolbox);

  /* get exercise */
  if (!toolbox_has_tool(exercise_toolbox, argv[parsed_args])) {
    error_set(err, "exercise '%s' not found; option -help lists possible "
                   "tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call exercise */
  if (!had_err) {
    if (!(toolfunc = toolbox_get(exercise_toolbox, argv[parsed_args]))) {
      tool = toolbox_get_tool(exercise_toolbox, argv[parsed_args]);
      assert(tool);
    }
    nargv = cstr_array_prefix_first(argv + parsed_args,
                                    error_get_progname(err));
    error_set_progname(err, nargv[0]);
    if (toolfunc)
      had_err = toolfunc(argc - parsed_args, (const char**) nargv, err);
    else
      had_err = tool_run(tool, argc - parsed_args, (const char**) nargv, err);
  }

  /* free */
  cstr_array_delete(nargv);

  return had_err;
}

Tool* gt_exercise(void)
{
  return tool_new(gt_exercise_arguments_new,
                  gt_exercise_arguments_delete,
                  gt_exercise_option_parser_new,
                  NULL,
                  gt_exercise_runner);
}
