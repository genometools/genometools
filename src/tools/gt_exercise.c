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

static void save_exercise_name(void *key, void *value, void *data)
{
  const char *exercisename;
  Array *exercisenames;
  assert(key && value && data);
  exercisename = (const char*) key;
  exercisenames = (Array*) data;
  array_add(exercisenames, exercisename);
}

static int show_exercise_tools(const char *progname, void *data, Error *err)
{
  Hashtable *exercise_tools;
  Array *exercisenames;
  unsigned long i;
  error_check(err);
  assert(data);
  exercise_tools = (Hashtable*) data;
  exercisenames = array_new(sizeof(const char*));
  hashtable_foreach(exercise_tools, save_exercise_name, exercisenames);
  printf("\nExercise tools:\n\n");
  assert(array_size(exercisenames));
  qsort(array_get_space(exercisenames), array_size(exercisenames),
        array_elem_size(exercisenames), compare);
  for (i = 0; i < array_size(exercisenames); i++) {
    puts(*(const char**) array_get(exercisenames, i));
  }
  array_free(exercisenames);
  return 0;
}

static OPrval parse_options(int *parsed_args, int argc, char **argv,
                            Hashtable *exercise_tools, Error *err)
{
  OptionParser *op;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] exercise_tool_name [argument ...]",
                         "Call exercise tool with name exercise_tool_name and "
                         "pass argument(s) to it.");
  option_parser_set_comment_func(op, show_exercise_tools, exercise_tools);
  oprval = option_parser_parse_min_args(op, parsed_args, argc, argv,
                                        versionfunc, 1, err);
  option_parser_free(op);
  return oprval;
}

void register_exercises(Hashtable *exercise_tools)
{
  assert(exercise_tools);
  hashtable_add(exercise_tools, "affinealign", gt_affinealign);
  hashtable_add(exercise_tools, "align", gt_align);
  hashtable_add(exercise_tools, "casino", gt_casino);
  hashtable_add(exercise_tools, "coin", gt_coin);
  hashtable_add(exercise_tools, "consensus_sa", gt_consensus_sa);
  hashtable_add(exercise_tools, "msaparse", gt_msaparse);
  hashtable_add(exercise_tools, "neighborjoining", gt_neighborjoining);
  hashtable_add(exercise_tools, "nussinov_rna_fold", gt_nussinov_rna_fold);
  hashtable_add(exercise_tools, "qgramdist", gt_qgramdist);
  hashtable_add(exercise_tools, "scorematrix", gt_scorematrix);
  hashtable_add(exercise_tools, "swalign", gt_swalign);
  hashtable_add(exercise_tools, "upgma", gt_upgma);
}

int gt_exercise(int argc, char *argv[], Error *err)
{
  Hashtable *exercise_tools;
  int (*exercise)(int, char**, Error*);
  int parsed_args, has_err = 0;
  char **nargv = NULL;
  error_check(err);

  /* option parsing */
  exercise_tools = hashtable_new(HASH_STRING, NULL, NULL);
  register_exercises(exercise_tools);
  switch (parse_options(&parsed_args, argc, argv, exercise_tools, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      hashtable_free(exercise_tools);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      hashtable_free(exercise_tools);
      return 0;
  }
  assert(parsed_args < argc);

  /* get exercise */
  if (!(exercise = hashtable_get(exercise_tools, argv[1]))) {
    error_set(err, "exercise '%s' not found; option -help lists possible tools",
              argv[1]);
    has_err = -1;
  }

  /* call exercise */
  if (!has_err) {
    nargv = cstr_array_prefix_first(argv+parsed_args, argv[0]);
    has_err = exercise(argc-parsed_args, nargv, err);
  }

  /* free */
  cstr_array_free(nargv);
  hashtable_free(exercise_tools);

  return has_err;
}
