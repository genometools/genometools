/*
  Copyright (c) 2003-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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
#include "gtt.h"
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "lfs.h"
#include "lpeg.h"
#include "md5.h"
#include "ldes56.h"
#include "libgtcore/cstr_array.h"
#include "libgtcore/ensure.h"
#include "libgtcore/fa.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/gtdatapath.h"
#include "libgtcore/log.h"
#include "libgtcore/ma.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xansi.h"
#include "libgtcore/yarandom.h"
#include "libgtext/gtdatahelp.h"
#include "libgtext/luahelper.h"
#include "libgtlua/gt_lua.h"
#include "libgtlua/interactive.h"

#ifdef LIBGTVIEW
#include "libgtview/luaconfig.h"
#endif

struct GTR {
  bool test,
       interactive,
       debug;
  unsigned int seed;
  Str *testspacepeak;
  Toolbox *tools;
  Hashtable *unit_tests;
  lua_State *L;
#ifdef LIBGTVIEW
  Config *config;
#endif
};

GTR* gtr_new(Error *err)
{
  GTR *gtr;
  int had_err = 0;
#ifdef LIBGTVIEW
  Str *config_file = NULL;
#endif
  gtr = ma_calloc(1, sizeof (GTR));
  gtr->testspacepeak = str_new();
  gtr->L = luaL_newstate();
  if (!gtr->L) {
    error_set(err, "out of memory (cannot create new lua state)");
    had_err = -1;
  }
  if (!had_err) {
    luaL_openlibs(gtr->L);    /* open the standard libraries */
    luaopen_gt(gtr->L);       /* open all GenomeTools libraries */
    luaopen_lfs(gtr->L);      /* open Lua filesystem */
    luaopen_lpeg(gtr->L);     /* open LPeg library */
    luaopen_md5_core(gtr->L); /* open MD5 library */
    luaopen_des56(gtr->L);
    had_err = lua_set_modules_path(gtr->L, err);
  }
#ifdef LIBGTVIEW
  if (!had_err) {
    if (!(gtr->config = config_new_with_state(gtr->L)))
      had_err = -1;
  }
  if (!had_err) {
    if (!(config_file = gtdata_get_path(error_get_progname(err), err)))
      had_err = -1;
  }
  if (!had_err) {
    str_append_cstr(config_file, "/config/view.lua");
    if (file_exists(str_get(config_file))) {
      if (config_load_file(gtr->config, config_file, err))
        had_err = -1;
      else
        lua_put_config_in_registry(gtr->L, gtr->config);
    }
  }
  str_delete(config_file);
#endif
  if (had_err) {
    ma_free(gtr);
    return NULL;
  }
  return gtr;
}

static int show_gtr_help(const char *progname, void *data, Error *err)
{
  int had_err;
  had_err = toolbox_show(progname, data, err);
  if (!had_err)
    had_err = gtdata_show_help(progname, NULL, err);
  return had_err;
}

OPrval gtr_parse(GTR *gtr, int *parsed_args, int argc, const char **argv,
                 Error *err)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;

  error_check(err);
  assert(gtr);
  op = option_parser_new("[option ...] [tool | script] [argument ...]",
                         "The GenomeTools (gt) genome analysis system "
                          "(http://genometools.org).");
  option_parser_set_comment_func(op, show_gtr_help, gtr->tools);
  o = option_new_bool("i", "enter interactive mode after executing 'tool' or "
                      "'script'", &gtr->interactive, false);
  option_hide_default(o);
  option_parser_add_option(op, o);
  o = option_new_bool("test", "perform unit tests and exit", &gtr->test, false);
  option_hide_default(o);
  option_parser_add_option(op, o);
  o = option_new_debug(&gtr->debug);
  option_parser_add_option(op, o);
  o = option_new_uint("seed", "set seed for random number generator manually\n"
                      "0 generates a seed from current time and process id",
                      &gtr->seed, 0);
  option_is_development_option(o);
  option_parser_add_option(op, o);
  o = option_new_filename("testspacepeak", "alloc 64 MB and mmap the given "
                          "file", gtr->testspacepeak);
  option_is_development_option(o);
  option_parser_add_option(op, o);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

void gtr_register_components(GTR *gtr)
{
  assert(gtr);
  /* add tools */
  toolbox_delete(gtr->tools);
  gtr->tools = gtt_tools();
  /* add unit tests */
  hashtable_delete(gtr->unit_tests);
  gtr->unit_tests = gtt_unit_tests();
}

int run_test(void *key, void *value, void *data, Error *e)
{
  const char *testname;
  int (*test)(Error*);
  int had_err, *had_errp;
  error_check(e);
  assert(key && value && data);
  testname = (const char*) key;
  test = (int (*)(Error*)) value;
  had_errp = (int*) data;
  printf("%s...", testname);
  xfflush(stdout);
  had_err = test(e);
  if (had_err) {
    xputs("error");
    *had_errp = had_err;
    fprintf(stderr, "first error: %s\n", error_get(e));
    error_unset(e);
    xfflush(stderr);
  }
  else
    xputs("ok");
  xfflush(stdout);
  return 0;
}

static int run_tests(GTR *gtr, Error *err)
{
  int test_err = 0, had_err = 0;
  error_check(err);
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

  /* show seed */
  printf("seed=%u\n", gtr->seed);

  if (gtr->unit_tests) {
    had_err = hashtable_foreach_ao(gtr->unit_tests, run_test, &test_err, err);
    assert(!had_err); /* cannot happen, run_test() is sane */
  }
  if (test_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

int gtr_run(GTR *gtr, int argc, const char **argv, Error *err)
{
  Toolfunc toolfunc;
  Tool *tool = NULL;
  char **nargv = NULL;
  void *mem, *map;
  int had_err = 0;
  error_check(err);
  assert(gtr);
  if (gtr->debug)
    log_enable();
  gtr->seed = ya_rand_init(gtr->seed);
  log_log("seed=%u", gtr->seed);
  if (gtr->test) {
    return run_tests(gtr, err);
  }
  if (str_length(gtr->testspacepeak)) {
    mem = ma_malloc(1 << 26); /* alloc 64 MB */;
    map = fa_mmap_read(str_get(gtr->testspacepeak), NULL);
    fa_xmunmap(map);
    ma_free(mem);
  }
  if (argc == 0 && !gtr->interactive) {
    error_set(err, "neither tool nor script specified; option -help lists "
                   "possible tools");
    had_err = -1;
  }
  if (!had_err && argc) {
    if (!gtr->tools || !toolbox_has_tool(gtr->tools, argv[0])) {
      /* no tool found -> try to open script */
      if (file_exists(argv[0])) {
        /* run script */
        nargv = cstr_array_prefix_first(argv, error_get_progname(err));
        lua_set_arg(gtr->L, nargv[0], (const char**) nargv+1);
        if (luaL_dofile(gtr->L, argv[0])) {
          /* error */
          assert(lua_isstring(gtr->L, -1)); /* error message on top */
          error_set(err, "could not execute script %s",
                    lua_tostring(gtr->L, -1));
          had_err = -1;
          lua_pop(gtr->L, 1); /* pop error message */
        }
      }
      else {
        /* neither tool nor script found */
        error_set(err, "neither tool nor script '%s' found; option -help lists "
                       "possible tools", argv[0]);
        had_err = -1;
      }
    }
    else {
      /* run tool */
      if (!(toolfunc = toolbox_get(gtr->tools, argv[0]))) {
        tool = toolbox_get_tool(gtr->tools, argv[0]);
        assert(tool);
      }
      nargv = cstr_array_prefix_first(argv, error_get_progname(err));
      error_set_progname(err, nargv[0]);
      if (toolfunc)
        had_err = toolfunc(argc, (const char**) nargv, err);
      else
        had_err = tool_run(tool, argc, (const char**) nargv, err);
    }
  }
  cstr_array_delete(nargv);
  if (!had_err && gtr->interactive) {
    showshortversion(error_get_progname(err));
    lua_set_arg(gtr->L, error_get_progname(err), argv);
    run_interactive_lua_interpreter(gtr->L);
  }
  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

void gtr_delete(GTR *gtr)
{
  if (!gtr) return;
  str_delete(gtr->testspacepeak);
  toolbox_delete(gtr->tools);
  hashtable_delete(gtr->unit_tests);
  if (gtr->L) lua_close(gtr->L);
#ifdef LIBGTVIEW
  config_delete_without_state(gtr->config);
#endif
  ma_free(gtr);
}
