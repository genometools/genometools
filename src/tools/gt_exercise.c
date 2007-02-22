/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"
#include "gt_affinealign.h"
#include "gt_align.h"
#include "gt_casino.h"
#include "gt_coin.h"
#include "gt_consensus_sa.h"
#include "gt_msaparse.h"
#include "gt_neighborjoining.h"
#include "gt_nussinov_rna_fold.h"
#include "gt_qgramdist.h"
#include "gt_scorematrix.h"
#include "gt_swalign.h"
#include "gt_upgma.h"

static int save_exercise_name(void *key, void *value, void *data, Env *env)
{
  const char *exercisename;
  Array *exercisenames;
  env_error_check(env);
  assert(key && value && data);
  exercisename = (const char*) key;
  exercisenames = (Array*) data;
  array_add(exercisenames, exercisename, env);
  return 0;
}

static int show_exercise_tools(const char *progname, void *data, Env *env)
{
  Hashtable *exercise_tools;
  Array *exercisenames;
  unsigned long i;
  int has_err;
  env_error_check(env);
  assert(data);
  exercise_tools = (Hashtable*) data;
  exercisenames = array_new(sizeof (const char*), env);
  has_err = hashtable_foreach(exercise_tools, save_exercise_name, exercisenames,
                              env);
  assert(!has_err); /* cannot happen, save_exercise_name() is sane */
  printf("\nExercise tools:\n\n");
  assert(array_size(exercisenames));
  qsort(array_get_space(exercisenames), array_size(exercisenames),
        array_elem_size(exercisenames), compare);
  for (i = 0; i < array_size(exercisenames); i++) {
    puts(*(const char**) array_get(exercisenames, i));
  }
  array_delete(exercisenames, env);
  return 0;
}

static OPrval parse_options(int *parsed_args, int argc, char **argv,
                            Hashtable *exercise_tools, Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] exercise_tool_name [argument ...]",
                         "Call exercise tool with name exercise_tool_name and "
                         "pass argument(s) to it.", env);
  option_parser_set_comment_func(op, show_exercise_tools, exercise_tools);
  oprval = option_parser_parse_min_args(op, parsed_args, argc, argv,
                                        versionfunc, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

void register_exercises(Hashtable *exercise_tools, Env *env)
{
  assert(exercise_tools);
  hashtable_add(exercise_tools, "affinealign", gt_affinealign, env);
  hashtable_add(exercise_tools, "align", gt_align, env);
  hashtable_add(exercise_tools, "casino", gt_casino, env);
  hashtable_add(exercise_tools, "coin", gt_coin, env);
  hashtable_add(exercise_tools, "consensus_sa", gt_consensus_sa, env);
  hashtable_add(exercise_tools, "msaparse", gt_msaparse, env);
  hashtable_add(exercise_tools, "neighborjoining", gt_neighborjoining, env);
  hashtable_add(exercise_tools, "nussinov_rna_fold", gt_nussinov_rna_fold, env);
  hashtable_add(exercise_tools, "qgramdist", gt_qgramdist, env);
  hashtable_add(exercise_tools, "scorematrix", gt_scorematrix, env);
  hashtable_add(exercise_tools, "swalign", gt_swalign, env);
  hashtable_add(exercise_tools, "upgma", gt_upgma, env);
}

int gt_exercise(int argc, char *argv[], Env *env)
{
  Hashtable *exercise_tools;
  int (*exercise)(int, char**, Env*);
  int parsed_args, has_err = 0;
  char **nargv = NULL;
  env_error_check(env);

  /* option parsing */
  exercise_tools = hashtable_new(HASH_STRING, NULL, NULL, env);
  register_exercises(exercise_tools, env);
  switch (parse_options(&parsed_args, argc, argv, exercise_tools, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      hashtable_delete(exercise_tools, env);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      hashtable_delete(exercise_tools, env);
      return 0;
  }
  assert(parsed_args < argc);

  /* get exercise */
  if (!(exercise = hashtable_get(exercise_tools, argv[1]))) {
    env_error_set(env, "exercise '%s' not found; option -help lists possible "
                  "tools", argv[1]);
    has_err = -1;
  }

  /* call exercise */
  if (!has_err) {
    nargv = cstr_array_prefix_first(argv+parsed_args, argv[0]);
    has_err = exercise(argc-parsed_args, nargv, env);
  }

  /* free */
  cstr_array_delete(nargv, env);
  hashtable_delete(exercise_tools, env);

  return has_err;
}
