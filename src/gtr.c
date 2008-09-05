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
#include "interactive.h"
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "lfs.h"
#include "lpeg.h"
#include "md5.h"
#include "ldes56.h"
#include "core/cstr_array.h"
#include "core/ensure.h"
#include "core/fa.h"
#include "core/fileutils.h"
#include "core/gtdatapath.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/versionfunc.h"
#include "core/xansi.h"
#include "core/yarandom.h"
#include "extended/feature_type_factory_builtin.h"
#include "extended/gtdatahelp.h"
#include "extended/luahelper.h"
#include "gtlua/gt_lua.h"

#ifndef WITHOUT_CAIRO
#include "annotationsketch/luastyle.h"
#include "annotationsketch/style_with_state.h"
#endif

struct GTR {
  bool test,
       interactive,
       debug,
       check64bit;
  unsigned int seed;
  Str *debugfp,
      *testspacepeak;
  Toolbox *tools;
  Hashmap *unit_tests;
  lua_State *L;
  FeatureTypeFactory *feature_type_factory; /* for gtlua */
#ifndef WITHOUT_CAIRO
  Style *style;
#endif
  FILE *logfp;
};

GTR* gtr_new(Error *err)
{
  GTR *gtr;
  int had_err = 0;
#ifndef WITHOUT_CAIRO
  Str *style_file = NULL;
#endif
  gtr = ma_calloc(1, sizeof (GTR));
  gtr->debugfp = str_new();
  gtr->testspacepeak = str_new();
  gtr->L = luaL_newstate();
  if (!gtr->L) {
    error_set(err, "out of memory (cannot create new lua state)");
    had_err = -1;
  }
  if (!had_err) {
    gtr->feature_type_factory = feature_type_factory_builtin_new();
    lua_put_feature_type_factory_in_registry(gtr->L, gtr->feature_type_factory);
    luaL_openlibs(gtr->L);    /* open the standard libraries */
    luaopen_gt(gtr->L);       /* open all GenomeTools libraries */
    luaopen_lfs(gtr->L);      /* open Lua filesystem */
    luaopen_lpeg(gtr->L);     /* open LPeg library */
    luaopen_md5_core(gtr->L); /* open MD5 library */
    luaopen_des56(gtr->L);
    had_err = lua_set_modules_path(gtr->L, err);
  }
#ifndef WITHOUT_CAIRO
  if (!had_err) {
    if (!(gtr->style = style_new_with_state(gtr->L)))
      had_err = -1;
  }
  if (!had_err) {
    if (!(style_file = gtdata_get_path(error_get_progname(err), err)))
      had_err = -1;
  }
  if (!had_err) {
    str_append_cstr(style_file, "/config/view.lua");
    if (file_exists(str_get(style_file))) {
      if (style_load_file(gtr->style, str_get(style_file), err))
        had_err = -1;
      else
        lua_put_style_in_registry(gtr->L, gtr->style);
    }
  }
  str_delete(style_file);
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
  Option *o, *debug_option, *debugfp_option;
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
  debug_option = option_new_debug(&gtr->debug);
  option_parser_add_option(op, debug_option);
  debugfp_option = option_new_string("debugfp",
                                     "set file pointer for debugging output\n"
                                     "use ``stdout'' for standard output\n"
                                     "use ``stderr'' for standard error\n"
                                     "or any other string to use the "
                                     "corresponding file (will be overwritten "
                                     "without warning!)", gtr->debugfp,
                                     "stderr");
  option_is_development_option(debugfp_option);
  option_parser_add_option(op, debugfp_option);
  option_imply(debugfp_option, debug_option);
  o = option_new_uint("seed", "set seed for random number generator manually\n"
                      "0 generates a seed from current time and process id",
                      &gtr->seed, 0);
  option_is_development_option(o);
  option_parser_add_option(op, o);
  o = option_new_bool("64bit", "exit with code 0 if this is a 64bit binary, "
                      "with 1 otherwise", &gtr->check64bit, false);
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
  hashmap_delete(gtr->unit_tests);
  gtr->unit_tests = gtt_unit_tests();
}

static int
run_test(void *key, void *value, void *data, Error *err)
{
  int had_err, *had_errp;
  char *testname = key;
  UnitTestFunc test = value;
  error_check(err);
  assert(testname && test && data);
  had_errp = (int*) data;
  printf("%s...", testname);
  xfflush(stdout);
  had_err = test(err);
  if (had_err) {
    xputs("error");
    *had_errp = had_err;
    fprintf(stderr, "first error: %s\n", error_get(err));
    error_unset(err);
    xfflush(stderr);
  }
  else
    xputs("ok");
  xfflush(stdout);
  return 0;
}

static int check64bit(void)
{
  if (sizeof (unsigned long) == 8)
    return EXIT_SUCCESS;
  return EXIT_FAILURE;
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
  hashmap_unit_test(err);
  if (gtr->unit_tests) {
    had_err = hashmap_foreach_in_key_order(
      gtr->unit_tests, run_test, &test_err, err);
    assert(!had_err); /* cannot happen, run_test() is sane */
  }
  if (test_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

static void enable_logging(const char *debugfp, FILE **logfp)
{
  log_enable();
  if (!strcmp(debugfp, "stdout"))
    log_set_fp(stdout);
  else if (!strcmp(debugfp, "stderr"))
    log_set_fp(stderr);
  else {
    *logfp = fa_xfopen(debugfp, "w");
    log_set_fp(*logfp);
  }
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
    enable_logging(str_get(gtr->debugfp), &gtr->logfp);
  gtr->seed = ya_rand_init(gtr->seed);
  log_log("seed=%u", gtr->seed);
  if (gtr->check64bit)
    return check64bit();
  if (gtr->test)
    return run_tests(gtr, err);
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
  fa_fclose(gtr->logfp);
  str_delete(gtr->testspacepeak);
  str_delete(gtr->debugfp);
  toolbox_delete(gtr->tools);
  hashmap_delete(gtr->unit_tests);
  feature_type_factory_delete(gtr->feature_type_factory);
  if (gtr->L) lua_close(gtr->L);
#ifndef WITHOUT_CAIRO
  style_delete_without_state(gtr->style);
#endif
  ma_free(gtr);
}
