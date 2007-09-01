/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/cstr.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/toolbox.h"
#include "gt_affinealign.h"
#include "gt_align.h"
#include "gt_casino.h"
#include "gt_coin.h"
#include "gt_consensus_sa.h"
#include "gt_linearalign.h"
#include "gt_msaparse.h"
#include "gt_multiset_matching.h"
#include "gt_neighborjoining.h"
#include "gt_nussinov_rna_fold.h"
#include "gt_qgramdist.h"
#include "gt_scorematrix.h"
#include "gt_swalign.h"
#include "gt_upgma.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Toolbox *exercise_toolbox, Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] exercise_tool_name [argument ...]",
                         "Call exercise tool with name exercise_tool_name and "
                         "pass argument(s) to it.", env);
  option_parser_set_comment_func(op, toolbox_show, exercise_toolbox);
  oprval = option_parser_parse_min_args(op, parsed_args, argc, argv,
                                        versionfunc, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

void register_exercises(Toolbox *exercise_toolbox, Env *env)
{
  assert(exercise_toolbox);
  toolbox_add(exercise_toolbox, "affinealign", gt_affinealign, env);
  toolbox_add(exercise_toolbox, "align", gt_align, env);
  toolbox_add(exercise_toolbox, "casino", gt_casino, env);
  toolbox_add(exercise_toolbox, "coin", gt_coin, env);
  toolbox_add(exercise_toolbox, "consensus_sa", gt_consensus_sa, env);
  toolbox_add(exercise_toolbox, "msaparse", gt_msaparse, env);
  toolbox_add(exercise_toolbox, "msmatch", gt_multiset_matching, env);
  toolbox_add(exercise_toolbox, "linearalign", gt_linearalign, env);
  toolbox_add(exercise_toolbox, "neighborjoining", gt_neighborjoining, env);
  toolbox_add(exercise_toolbox, "nussinov_rna_fold", gt_nussinov_rna_fold, env);
  toolbox_add(exercise_toolbox, "qgramdist", gt_qgramdist, env);
  toolbox_add(exercise_toolbox, "scorematrix", gt_scorematrix, env);
  toolbox_add(exercise_toolbox, "swalign", gt_swalign, env);
  toolbox_add(exercise_toolbox, "upgma", gt_upgma, env);
}

int gt_exercise(int argc, const char **argv, Env *env)
{
  Toolbox *exercise_toolbox;
  Tool exercise;
  int parsed_args, had_err = 0;
  char **nargv = NULL;
  env_error_check(env);

  /* option parsing */
  exercise_toolbox = toolbox_new(env);
  register_exercises(exercise_toolbox, env);
  switch (parse_options(&parsed_args, argc, argv, exercise_toolbox, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      toolbox_delete(exercise_toolbox, env);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      toolbox_delete(exercise_toolbox, env);
      return 0;
  }
  assert(parsed_args < argc);

  /* get exercise */
  if (!(exercise = toolbox_get(exercise_toolbox, argv[1]))) {
    env_error_set(env, "exercise '%s' not found; option -help lists possible "
                  "tools", argv[1]);
    had_err = -1;
  }

  /* call exercise */
  if (!had_err) {
    nargv = cstr_array_prefix_first(argv+parsed_args, argv[0], env);
    had_err = exercise(argc-parsed_args, (const char**) nargv, env);
  }

  /* free */
  cstr_array_delete(nargv, env);
  toolbox_delete(exercise_toolbox, env);

  return had_err;
}
