/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"
#include "gtr.h"
#include "lua.h"
#include "lauxlib.h"
#include "tools/gt_bioseq.h"
#include "tools/gt_cds.h"
#include "tools/gt_clean.h"
#include "tools/gt_csa.h"
#include "tools/gt_eval.h"
#include "tools/gt_exercise.h"
#include "tools/gt_extractfeat.h"
#include "tools/gt_filter.h"
#include "tools/gt_gff3.h"
#include "tools/gt_gtf2gff3.h"
#include "tools/gt_merge.h"
#include "tools/gt_mmapandread.h"
#include "tools/gt_mutate.h"
#include "tools/gt_stat.h"

struct GTR {
  bool test,
       interactive,
       debug;
  Hashtable *tools,
            *unit_tests;
  lua_State *L;
};

GTR* gtr_new(Env *env)
{
  GTR *gtr = env_ma_calloc(env, 1, sizeof (GTR));
  /* XXX */
  gtr->L = luaL_newstate();
  assert(gtr->L);
  return gtr;
}

static int show_tool(void *key, void *value, void *data, Env *env)
{
  const char *toolname;
  Array *toolnames;
  env_error_check(env);
  assert(key && value && data);
  toolname = (const char*) key;
  toolnames = (Array*) data;
  array_add(toolnames, toolname, env);
  return 0;
}

static int show_option_comments(const char *progname, void *data, Env *env)
{
  Array *toolnames;
  unsigned long i;
  int has_err;
  GTR *gtr;
  env_error_check(env);
  assert(data);
  gtr = (GTR*) data;
  toolnames = array_new(sizeof (const char*), env);
  if (gtr->tools) {
    has_err = hashtable_foreach(gtr->tools, show_tool, toolnames, env);
    assert(!has_err); /* cannot happen, show_tool() is sane */
    printf("\nTools:\n\n");
    assert(array_size(toolnames));
    qsort(array_get_space(toolnames), array_size(toolnames),
          array_elem_size(toolnames), compare);
    for (i = 0; i < array_size(toolnames); i++)
      xputs(*(const char**) array_get(toolnames, i));
  }
  array_delete(toolnames, env);
  return 0;
}

OPrval gtr_parse(GTR *gtr, int *parsed_args, int argc, const char **argv,
                 Env *env)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  env_error_check(env);
  assert(gtr);
  op = option_parser_new("[option ...] [tool ...] [argument ...]",
                         "The GenomeTools (gt) genome analysis system "
                          "(http://genometools.org).", env);
  option_parser_set_comment_func(op, show_option_comments, gtr);
  o = option_new_bool("test", "perform unit tests and exit", &gtr->test, false,
                      env);
  option_parser_add_option(op, o, env);
  o = option_new_bool("i", "enter interactive mode after executing 'tool'",
                      &gtr->interactive, false, env);
  option_is_development_option(o);
  option_parser_add_option(op, o, env);
  o = option_new_debug(&gtr->debug, env);
  option_parser_add_option(op, o, env);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  option_parser_delete(op, env);
  (*parsed_args)--;
  return oprval;
}

void gtr_register_components(GTR *gtr, Env *env)
{
  assert(gtr);
  /* add tools */
  hashtable_delete(gtr->tools, env);
  gtr->tools = hashtable_new(HASH_STRING, NULL, NULL, env);
  hashtable_add(gtr->tools, "bioseq", gt_bioseq, env);
  hashtable_add(gtr->tools, "cds", gt_cds, env);
  hashtable_add(gtr->tools, "clean", gt_clean, env);
  hashtable_add(gtr->tools, "csa", gt_csa, env);
  hashtable_add(gtr->tools, "eval", gt_eval, env);
  hashtable_add(gtr->tools, "exercise", gt_exercise, env);
  hashtable_add(gtr->tools, "extractfeat", gt_extractfeat, env);
  hashtable_add(gtr->tools, "filter", gt_filter, env);
  hashtable_add(gtr->tools, "gff3", gt_gff3, env);
  hashtable_add(gtr->tools, "gtf2gff3", gt_gtf2gff3, env);
  hashtable_add(gtr->tools, "merge", gt_merge, env);
  hashtable_add(gtr->tools, "mmapandread", gt_mmapandread, env);
  hashtable_add(gtr->tools, "mutate", gt_mutate, env);
  hashtable_add(gtr->tools, "stat", gt_stat, env);
  /* add unit tests */
  hashtable_delete(gtr->unit_tests, env);
  gtr->unit_tests = hashtable_new(HASH_STRING, NULL, NULL, env);
  hashtable_add(gtr->unit_tests, "alignment class", alignment_unit_test, env);
  hashtable_add(gtr->unit_tests, "array class", array_unit_test, env);
  hashtable_add(gtr->unit_tests, "bittab class", bittab_unit_test, env);
  hashtable_add(gtr->unit_tests, "bsearch module", bsearch_unit_test, env);
  hashtable_add(gtr->unit_tests, "countingsort module", countingsort_unit_test,
                env);
  hashtable_add(gtr->unit_tests, "dlist class", dlist_unit_test, env);
  hashtable_add(gtr->unit_tests, "evaluator class", evaluator_unit_test, env);
  hashtable_add(gtr->unit_tests, "grep module", grep_unit_test, env);
  hashtable_add(gtr->unit_tests, "hashtable class", hashtable_unit_test, env);
  hashtable_add(gtr->unit_tests, "hmm class", hmm_unit_test, env);
  hashtable_add(gtr->unit_tests, "range class", range_unit_test, env);
  hashtable_add(gtr->unit_tests, "splicedseq class", splicedseq_unit_test, env);
  hashtable_add(gtr->unit_tests, "splitter class", splitter_unit_test, env);
  hashtable_add(gtr->unit_tests, "string class", str_unit_test, env);
  hashtable_add(gtr->unit_tests, "tokenizer class", tokenizer_unit_test, env);
}

int run_test(void *key, void *value, void *data, Env *env)
{
  const char *testname;
  int (*test)(Env*);
  int has_err, *had_errp;
  env_error_check(env);
  assert(key && value && data);
  testname = (const char*) key;
  test = value;
  had_errp = (int*) data;
  printf("%s...", testname);
  xfflush(stdout);
  has_err = test(env);
  if (has_err) {
    xputs("error");
    *had_errp = has_err;
    fprintf(stderr, "first error: %s\n", env_error_get(env));
    env_error_unset(env);
    xfflush(stderr);
  }
  else
    xputs("ok");
  xfflush(stdout);
  return 0;
}

static int run_tests(GTR *gtr, Env *env)
{
  int had_err = 0, has_err = 0;
  env_error_check(env);
  assert(gtr);

  /* The following type assumptions are made in the GenomeTools library. */
  ensure(has_err, sizeof (char) == 1);
  ensure(has_err, sizeof (unsigned char) == 1);
  ensure(has_err, sizeof (short) == 2);
  ensure(has_err, sizeof (unsigned short) == 2);
  ensure(has_err, sizeof (int) == 4);
  ensure(has_err, sizeof (unsigned int) == 4);
  ensure(has_err, sizeof (long) == 4 || sizeof (long) == 8);
  ensure(has_err, sizeof (unsigned long) == 4 || sizeof (unsigned long) == 8);
  ensure(has_err, sizeof (unsigned long) >= sizeof (size_t));
  ensure(has_err, sizeof (long long) == 8);
  ensure(has_err, sizeof (unsigned long long) == 8);

  if (gtr->unit_tests) {
    has_err = hashtable_foreach(gtr->unit_tests, run_test, &had_err, env);
    assert(!has_err); /* cannot happen, run_test() is sane */
  }
  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

int gtr_run(GTR *gtr, int argc, const char **argv, Env *env)
{
  int (*tool)(int, char**, Env*) = NULL;
  char **nargv = NULL;
  int has_err = 0;
  env_error_check(env);
  assert(gtr);
  if (gtr->test) {
    return run_tests(gtr, env);
  }
  if (gtr->debug)
    env_set_log(env, log_new(env_ma(env)));
  assert(argc);
  if (argc == 1 && !gtr->interactive) {
    env_error_set(env, "no tool specified; option -help lists possible tools");
    has_err = -1;
  }
  if (!has_err && argc > 1) {
    if (!gtr->tools || !(tool = hashtable_get(gtr->tools, argv[1]))) {
      env_error_set(env, "tool '%s' not found; option -help lists possible "
                         "tools", argv[1]);
      has_err = -1;
    }
  }
  if (!has_err && argc > 1) {
    nargv = cstr_array_prefix_first(argv+1, argv[0], env);
    has_err = tool(argc-1, nargv, env);
  }
  cstr_array_delete(nargv, env);
  if (!has_err && gtr->interactive) {
    /* XXX */
    env_error_set(env, "interactive mode not implemented yet");
    has_err = -1;
  }
  if (has_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

void gtr_delete(GTR *gtr, Env *env)
{
  if (!gtr) return;
  hashtable_delete(gtr->tools, env);
  hashtable_delete(gtr->unit_tests, env);
  if (gtr->L) lua_close(gtr->L);
  env_ma_free(gtr, env);
}
