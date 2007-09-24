/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg

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

#include "gtr.h"
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "lfs.h"
#include "libgtcore/array.h"
#include "libgtcore/bitpackarray.h"
#include "libgtcore/bitpackstring.h"
#include "libgtcore/bittab.h"
#include "libgtcore/countingsort.h"
#include "libgtcore/cstr.h"
#include "libgtcore/dlist.h"
#include "libgtcore/dynbittab.h"
#include "libgtcore/ensure.h"
#include "libgtcore/env.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/getbasename.h"
#include "libgtcore/gtdatapath.h"
#include "libgtcore/grep.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/range.h"
#include "libgtcore/splitter.h"
#include "libgtcore/tokenizer.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xansi.h"
#include "libgtext/alignment.h"
#include "libgtext/bsearch.h"
#include "libgtext/evaluator.h"
#include "libgtext/gtdatahelp.h"
#include "libgtext/hmm.h"
#include "libgtext/splicedseq.h"
#include "libgtext/toolbox.h"
#include "libgtlua/gt_lua.h"
#include "libgtlua/helper.h"
#include "libgtlua/interactive.h"
#include "tools/gt_bioseq.h"
#include "tools/gt_cds.h"
#include "tools/gt_chseqids.h"
#include "tools/gt_clean.h"
#include "tools/gt_csa.h"
#include "tools/gt_dev.h"
#include "tools/gt_eval.h"
#include "tools/gt_exercise.h"
#include "tools/gt_extractfeat.h"
#include "tools/gt_filter.h"
#include "tools/gt_gff3.h"
#include "tools/gt_gff3_to_gtf.h"
#include "tools/gt_gtf_to_gff3.h"
#include "tools/gt_ltrharvest.h"
#include "tools/gt_merge.h"
#include "tools/gt_mmapandread.h"
#include "tools/gt_mutate.h"
#include "tools/gt_splitfasta.h"
#include "tools/gt_stat.h"
#include "tools/gt_suffixerator.h"
#include "tools/gt_mkfmindex.h"
#include "tools/gt_uniquesub.h"

#ifdef LIBGTVIEW
#include "libgtview/block.h"
#include "libgtview/config.h"
#include "libgtview/diagram.h"
#include "libgtview/feature_index.h"
#include "libgtview/gt_view.h"
#include "libgtview/track.h"
#endif

struct GTR {
  bool test,
       interactive,
       debug;
  Str *testspacepeak;
  Toolbox *toolbox;
  Hashtable *unit_tests;
  lua_State *L;
#ifdef LIBGTVIEW
  Config *config;
#endif
};

GTR* gtr_new(Env *env)
{
  GTR *gtr = env_ma_calloc(env, 1, sizeof (GTR));
#ifdef LIBGTVIEW
  Str *config_file;
#endif
  gtr->testspacepeak = str_new(env);
  gtr->L = luaL_newstate();
  assert(gtr->L); /* XXX: proper error message  */
  luaL_openlibs(gtr->L); /* open the standard libraries */
  put_env_in_registry(gtr->L, env); /* we have to register the env object,
                                       before we can open the GenomeTools
                                       libraries */
  luaopen_gt(gtr->L); /* open all GenomeTools libraries */
  luaopen_lfs(gtr->L); /* open Lua filesystem */
#ifdef LIBGTVIEW
  gtr->config = config_new_with_state(gtr->L, env);
  config_file = gtdata_get_path(env_error_get_progname(env), env);
  str_append_cstr(config_file, "/config/view.lua", env);
  if (file_exists(str_get(config_file))) {
    if (config_load_file(gtr->config, config_file, env)) {
      /* XXX: hack... */
      fprintf(stderr, "%s: error: %s\n", env_error_get_progname(env),
              env_error_get(env));
      exit(EXIT_FAILURE);
    }
    else
      put_config_in_registry(gtr->L, gtr->config);
  }
  str_delete(config_file, env);
#endif
  return gtr;
}

static int show_gtr_help(const char *progname, void *data, Env* env)
{
  int had_err;
  had_err = toolbox_show(progname, data, env);
  if (!had_err)
    had_err = gtdata_show_help(progname, NULL, env);
  return had_err;
}

OPrval gtr_parse(GTR *gtr, int *parsed_args, int argc, const char **argv,
                 Env *env)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;

  env_error_check(env);
  assert(gtr);
  op = option_parser_new("[option ...] [tool | script] [argument ...]",
                         "The GenomeTools (gt) genome analysis system "
                          "(http://genometools.org).", env);
  option_parser_set_comment_func(op, show_gtr_help, gtr->toolbox);
  o = option_new_bool("i", "enter interactive mode after executing 'tool' or "
                      "'script'", &gtr->interactive, false, env);
  option_hide_default(o);
  option_parser_add_option(op, o, env);
  o = option_new_bool("test", "perform unit tests and exit", &gtr->test, false,
                      env);
  option_hide_default(o);
  option_parser_add_option(op, o, env);
  o = option_new_debug(&gtr->debug, env);
  option_parser_add_option(op, o, env);
  o = option_new_filename("testspacepeak", "alloc 64 MB and mmap the given "
                          "file", gtr->testspacepeak, env);
  option_is_development_option(o);
  option_parser_add_option(op, o, env);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  option_parser_delete(op, env);
  return oprval;
}

void gtr_register_components(GTR *gtr, Env *env)
{
  assert(gtr);
  /* add tools */
  toolbox_delete(gtr->toolbox, env);
  gtr->toolbox = toolbox_new(env);
  toolbox_add(gtr->toolbox, "bioseq", gt_bioseq, env);
  toolbox_add(gtr->toolbox, "cds", gt_cds, env);
  toolbox_add(gtr->toolbox, "chseqids", gt_chseqids, env);
  toolbox_add(gtr->toolbox, "clean", gt_clean, env);
  toolbox_add(gtr->toolbox, "csa", gt_csa, env);
  toolbox_add(gtr->toolbox, "dev", gt_dev, env);
  toolbox_add(gtr->toolbox, "eval", gt_eval, env);
  toolbox_add(gtr->toolbox, "exercise", gt_exercise, env);
  toolbox_add(gtr->toolbox, "extractfeat", gt_extractfeat, env);
  toolbox_add(gtr->toolbox, "filter", gt_filter, env);
  toolbox_add(gtr->toolbox, "gff3", gt_gff3, env);
  toolbox_add(gtr->toolbox, "gff3_to_gtf", gt_gff3_to_gtf, env);
  toolbox_add(gtr->toolbox, "gtf_to_gff3", gt_gtf_to_gff3, env);
  toolbox_add(gtr->toolbox, "ltrharvest", gt_ltrharvest, env);
  toolbox_add(gtr->toolbox, "merge", gt_merge, env);
  toolbox_add(gtr->toolbox, "mmapandread", gt_mmapandread, env);
  toolbox_add(gtr->toolbox, "mutate", gt_mutate, env);
  toolbox_add(gtr->toolbox, "splitfasta", gt_splitfasta, env);
  toolbox_add(gtr->toolbox, "stat", gt_stat, env);
  toolbox_add(gtr->toolbox, "suffixerator", gt_suffixerator, env);
  toolbox_add(gtr->toolbox, "mkfmindex", gt_mkfmindex, env);
  toolbox_add(gtr->toolbox, "uniquesub", gt_uniquesub, env);
#ifdef LIBGTVIEW
  toolbox_add(gtr->toolbox, "view", gt_view, env);
#endif
  /* add unit tests */
  hashtable_delete(gtr->unit_tests, env);
  gtr->unit_tests = hashtable_new(HASH_STRING, NULL, NULL, env);
  hashtable_add(gtr->unit_tests, "alignment class", alignment_unit_test, env);
  hashtable_add(gtr->unit_tests, "array class", array_unit_test, env);
  hashtable_add(gtr->unit_tests, "array example", array_example, env);
  hashtable_add(gtr->unit_tests, "bit pack array class",
                bitPackArray_unit_test, env);
  hashtable_add(gtr->unit_tests, "bit pack string module",
                bitPackString_unit_test, env);
  hashtable_add(gtr->unit_tests, "bittab class", bittab_unit_test, env);
  hashtable_add(gtr->unit_tests, "bsearch module", bsearch_unit_test, env);
  hashtable_add(gtr->unit_tests, "countingsort module", countingsort_unit_test,
                env);
  hashtable_add(gtr->unit_tests, "dlist class", dlist_unit_test, env);
  hashtable_add(gtr->unit_tests, "dynamic bittab class", dynbittab_unit_test,
                env);
  hashtable_add(gtr->unit_tests, "evaluator class", evaluator_unit_test, env);
  hashtable_add(gtr->unit_tests, "getbasename module", getbasename_unit_test,
                env);
  hashtable_add(gtr->unit_tests, "grep module", grep_unit_test, env);
  hashtable_add(gtr->unit_tests, "hashtable class", hashtable_unit_test, env);
  hashtable_add(gtr->unit_tests, "hmm class", hmm_unit_test, env);
  hashtable_add(gtr->unit_tests, "range class", range_unit_test, env);
  hashtable_add(gtr->unit_tests, "splicedseq class", splicedseq_unit_test, env);
  hashtable_add(gtr->unit_tests, "splitter class", splitter_unit_test, env);
  hashtable_add(gtr->unit_tests, "string class", str_unit_test, env);
  hashtable_add(gtr->unit_tests, "tokenizer class", tokenizer_unit_test, env);
#ifdef LIBGTVIEW
  hashtable_add(gtr->unit_tests, "block class", block_unit_test, env);
  hashtable_add(gtr->unit_tests, "config class", config_unit_test, env);
  hashtable_add(gtr->unit_tests, "diagram class", diagram_unit_test, env);
  hashtable_add(gtr->unit_tests, "element class", element_unit_test, env);
  hashtable_add(gtr->unit_tests, "feature index class", feature_index_unit_test,
                env);
  hashtable_add(gtr->unit_tests, "line class", line_unit_test, env);
  hashtable_add(gtr->unit_tests, "track class", track_unit_test, env);
#endif
}

int run_test(void *key, void *value, void *data, Env *env)
{
  const char *testname;
  int (*test)(Env*);
  int had_err, *had_errp;
  env_error_check(env);
  assert(key && value && data);
  testname = (const char*) key;
  test = (int (*)(Env *)) value;
  had_errp = (int*) data;
  printf("%s...", testname);
  xfflush(stdout);
  had_err = test(env);
  if (had_err) {
    xputs("error");
    *had_errp = had_err;
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
  int test_err = 0, had_err = 0;
  env_error_check(env);
  assert(gtr);

  /* The following type assumptions are made in the GenomeTools library. */
  ensure(had_err, sizeof (char) == 1);
  ensure(had_err, sizeof (unsigned char) == 1);
  ensure(had_err, sizeof (short) == 2);
  ensure(had_err, sizeof (unsigned short) == 2);
  ensure(had_err, sizeof (int) == 4);
  ensure(had_err, sizeof (unsigned int) == 4);
  ensure(had_err, sizeof (long) == 4 || sizeof (long) == 8);
  ensure(had_err, sizeof (unsigned long) == 4 || sizeof (unsigned long) == 8);
  ensure(had_err, sizeof (unsigned long) >= sizeof (size_t));
  ensure(had_err, sizeof (long long) == 8);
  ensure(had_err, sizeof (unsigned long long) == 8);

  if (gtr->unit_tests) {
    had_err = hashtable_foreach_ao(gtr->unit_tests, run_test, &test_err, env);
    assert(!had_err); /* cannot happen, run_test() is sane */
  }
  if (test_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

int gtr_run(GTR *gtr, int argc, const char **argv, Env *env)
{
  Tool tool = NULL;
  char **nargv = NULL;
  void *mem, *map;
  int had_err = 0;
  env_error_check(env);
  assert(gtr);
  if (gtr->debug)
    env_set_log(env, log_new(env_ma(env)));
  if (gtr->test) {
    return run_tests(gtr, env);
  }
  if (str_length(gtr->testspacepeak)) {
    mem = env_ma_malloc(env, 1 << 26); /* alloc 64 MB */;
    map = env_fa_mmap_read(env, str_get(gtr->testspacepeak), NULL);
    env_fa_xmunmap(map, env);
    env_ma_free(mem, env);
  }
  if (argc == 0 && !gtr->interactive) {
    env_error_set(env, "neither tool nor script specified; option -help lists "
                       "possible tools");
    had_err = -1;
  }
  if (!had_err && argc) {
    if (!gtr->toolbox || !(tool = toolbox_get(gtr->toolbox, argv[0]))) {
      /* no tool found -> try to open script */
      if (file_exists(argv[0])) {
        /* run script */
        nargv = cstr_array_prefix_first(argv, env_error_get_progname(env), env);
        set_arg_in_lua_interpreter(gtr->L, nargv[0], (const char**) nargv+1);
        if (luaL_dofile(gtr->L, argv[0])) {
          /* error */
          assert(lua_isstring(gtr->L, -1)); /* error message on top */
          env_error_set(env, "could not execute script %s",
                        lua_tostring(gtr->L, -1));
          had_err = -1;
          lua_pop(gtr->L, 1); /* pop error message */
        }
      }
      else {
        /* neither tool nor script found */
        env_error_set(env, "neither tool nor script '%s' found; option -help "
                           "lists possible tools", argv[0]);
        had_err = -1;
      }
    }
    else {
      /* run tool */
      nargv = cstr_array_prefix_first(argv, env_error_get_progname(env), env);
      env_error_set_progname(env, nargv[0]);
      had_err = tool(argc, (const char**) nargv, env);
    }
  }
  cstr_array_delete(nargv, env);
  if (!had_err && gtr->interactive) {
    showshortversion(env_error_get_progname(env));
    set_arg_in_lua_interpreter(gtr->L, env_error_get_progname(env), argv);
    run_interactive_lua_interpreter(gtr->L, env);
  }
  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

void gtr_delete(GTR *gtr, Env *env)
{
  if (!gtr) return;
  str_delete(gtr->testspacepeak, env);
  toolbox_delete(gtr->toolbox, env);
  hashtable_delete(gtr->unit_tests, env);
  if (gtr->L) lua_close(gtr->L);
#ifdef LIBGTVIEW
  config_delete_without_state(gtr->config, env);
#endif
  env_ma_free(gtr, env);
}
