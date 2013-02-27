/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <string.h>
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
#include "core/fileutils_api.h"
#include "core/gtdatapath.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/thread_api.h"
#include "core/unit_testing.h"
#include "core/versionfunc.h"
#include "core/xansi_api.h"
#include "core/yarandom.h"
#include "extended/gtdatahelp.h"
#include "extended/luahelper.h"
#include "extended/toolbox.h"
#include "gtlua/gt_lua.h"

#ifndef WITHOUT_CAIRO
#include "annotationsketch/luastyle.h"
#include "annotationsketch/style.h"
#endif

struct GtR {
  bool test,
       interactive,
       debug,
       check64bit;
  unsigned int seed;
  GtStr *debugfp,
      *testspacepeak;
  GtToolbox *tools;
  GtHashmap *unit_tests;
  GtStr *test_only;
  lua_State *L;
#ifndef WITHOUT_CAIRO
  GtStyle *style;
#endif
  FILE *logfp;
};

GtR* gtr_new(GtError *err)
{
  GtR *gtr;
  int had_err = 0;
#ifndef WITHOUT_CAIRO
  GtStr *style_file = NULL;
#endif
  gtr = gt_calloc(1, sizeof (GtR));
  gtr->debugfp = gt_str_new();
  gtr->testspacepeak = gt_str_new();
  gtr->test_only = gt_str_new();
  gtr->L = luaL_newstate();
  if (!gtr->L) {
    gt_error_set(err, "out of memory (cannot create new lua state)");
    had_err = -1;
  }
  if (!had_err) {
    luaL_openlibs(gtr->L);    /* open the standard libraries */
    gt_lua_open_lib(gtr->L);  /* open the GenomeTools library */
    lua_pushcfunction(gtr->L, luaopen_lpeg);
    lua_pushstring(gtr->L, "lpeg");
    lua_call(gtr->L, 1, 0);   /* open LPeg library */
    lua_pushcfunction(gtr->L, luaopen_md5_core);
    lua_pushstring(gtr->L, "md5");
    lua_call(gtr->L, 1, 0);   /* open MD5 library */
    lua_pushcfunction(gtr->L, luaopen_lfs);
    lua_pushstring(gtr->L, "lfs");
    lua_call(gtr->L, 1, 0);   /* open Lua filesystem */
    lua_pushcfunction(gtr->L, luaopen_des56);
    lua_pushstring(gtr->L, "des56");
    lua_call(gtr->L, 1, 0);   /* open DES56 library */
    had_err = gt_lua_set_modules_path(gtr->L, err);
  }
#ifndef WITHOUT_CAIRO
  if (!had_err) {
    lua_settop(gtr->L, 0);
    if (!(gtr->style = gt_style_new_with_state(gtr->L)))
      had_err = -1;
  }
  if (!had_err) {
    if (!(style_file = gt_get_gtdata_path(gt_error_get_progname(err), err)))
      had_err = -1;
  }
  if (!had_err) {
    gt_str_append_cstr(style_file, "/sketch/default.style");
    if (gt_file_exists(gt_str_get(style_file))) {
      if (gt_style_load_file(gtr->style, gt_str_get(style_file), err))
        had_err = -1;
      else
        gt_lua_put_style_in_registry(gtr->L, gtr->style);
    }
  }
  gt_str_delete(style_file);
#endif
  if (had_err) {
    gt_free(gtr);
    return NULL;
  }
  return gtr;
}

static int show_gtr_help(const char *progname, void *data, GtError *err)
{
  int had_err;
  had_err = gt_toolbox_show(progname, data, err);
  if (!had_err)
    had_err = gt_gtdata_show_help(progname, NULL, err);
  return had_err;
}

GtOPrval gtr_parse(GtR *gtr, int *parsed_args, int argc, const char **argv,
                   GtError *err)
{
  GtOptionParser *op;
  GtOption *o, *only_option, *debug_option, *debugfp_option;
  GtOPrval oprval;

  gt_error_check(err);
  gt_assert(gtr);
  op = gt_option_parser_new("[option ...] [tool | script] [argument ...]",
                            "The GenomeTools (gt) genome analysis system "
                            "(http://genometools.org).");
  gt_option_parser_set_comment_func(op, show_gtr_help, gtr->tools);
  o = gt_option_new_bool("i",
                         "enter interactive mode after executing 'tool' or "
                         "'script'", &gtr->interactive, false);
  gt_option_hide_default(o);
  gt_option_parser_add_option(op, o);
  o = gt_option_new_uint_min("j", "set number of parallel threads used at once",
                             &gt_jobs, 1, 1);
  gt_option_is_development_option(o);
  gt_option_parser_add_option(op, o);
  o = gt_option_new_bool("test", "perform unit tests and exit", &gtr->test,
                         false);
  gt_option_hide_default(o);
  gt_option_parser_add_option(op, o);
  only_option = gt_option_new_string("only", "perform single unit test "
                                     "(requires -test)", gtr->test_only, "");
  gt_option_imply(only_option, o);
  gt_option_is_development_option(only_option);
  gt_option_hide_default(only_option);
  gt_option_parser_add_option(op, only_option);
  debug_option = gt_option_new_debug(&gtr->debug);
  gt_option_parser_add_option(op, debug_option);
  debugfp_option = gt_option_new_string("debugfp",
                                     "set file pointer for debugging output\n"
                                     "use ``stdout'' for standard output\n"
                                     "use ``stderr'' for standard error\n"
                                     "or any other string to use the "
                                     "corresponding file (will be overwritten "
                                     "without warning!)", gtr->debugfp,
                                     "stderr");
  gt_option_is_development_option(debugfp_option);
  gt_option_parser_add_option(op, debugfp_option);
  gt_option_imply(debugfp_option, debug_option);
  o = gt_option_new_uint("seed",
                         "set seed for random number generator manually\n"
                         "0 generates a seed from current time and process id",
                         &gtr->seed, 0);
  gt_option_is_development_option(o);
  gt_option_parser_add_option(op, o);
  o = gt_option_new_bool("64bit", "exit with code 0 if this is a 64bit binary, "
                         "with 1 otherwise", &gtr->check64bit, false);
  gt_option_is_development_option(o);
  gt_option_parser_add_option(op, o);
  o = gt_option_new_filename("testspacepeak", "alloc 64 MB and mmap the given "
                             "file", gtr->testspacepeak);
  gt_option_is_development_option(o);
  gt_option_parser_add_option(op, o);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

void gtr_register_components(GtR *gtr)
{
  gt_assert(gtr);
  /* add tools */
  gt_toolbox_delete(gtr->tools);
  gtr->tools = gtt_tools();
  /* add unit tests */
  gt_hashmap_delete(gtr->unit_tests);
  gtr->unit_tests = gtt_unit_tests();
}

static int check64bit(void)
{
  if (sizeof (unsigned long) == 8)
    return EXIT_SUCCESS;
  return EXIT_FAILURE;
}

static int run_tests(GtR *gtr, GtError *err)
{
  int test_err = 0, had_err = 0;
  char* key;
  void* value;
  gt_error_check(err);
  gt_assert(gtr);

  /* The following type assumptions are made in the GenomeTools library. */
  gt_ensure(had_err, sizeof (char) == 1);
  gt_ensure(had_err, sizeof (unsigned char) == 1);
  gt_ensure(had_err, sizeof (short) == 2);
  gt_ensure(had_err, sizeof (unsigned short) == 2);
  gt_ensure(had_err, sizeof (int) == 4);
  gt_ensure(had_err, sizeof (unsigned int) == 4);
  gt_ensure(had_err, sizeof (long) == 4 || sizeof (long) == 8);
  gt_ensure(had_err, sizeof (unsigned long) == 4 ||
                     sizeof (unsigned long) == 8);
  gt_ensure(had_err, sizeof (unsigned long) >= sizeof (size_t));
  gt_ensure(had_err, sizeof (long long) == 8);
  gt_ensure(had_err, sizeof (unsigned long long) == 8);

  /* show seed */
  printf("seed=%u\n", gtr->seed);
  gt_hashmap_unit_test(err);
  if (gtr->unit_tests) {
    if (gt_str_length(gtr->test_only) > 0) {
      key = gt_str_get(gtr->test_only);
      value = gt_hashmap_get(gtr->unit_tests, key);
      if (value) {
        had_err = gt_unit_test_run(key, value, &test_err, err);
        gt_assert(!had_err); /* cannot happen, gt_unit_test_run() is sane */
      }
      else {
        gt_error_set(err, "Test \"%s\" not found", key);
        return EXIT_FAILURE;
      }
    }
    else {
      had_err = gt_hashmap_foreach_in_key_order(gtr->unit_tests,
                                                gt_unit_test_run, &test_err,
                                                err);
      gt_assert(!had_err); /* cannot happen, gt_unit_test_run() is sane */
    }
  }
  if (test_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

static void enable_logging(const char *debugfp, FILE **logfp)
{
  gt_log_enable();
  if (!strcmp(debugfp, "stdout"))
    gt_log_set_fp(stdout);
  else if (!strcmp(debugfp, "stderr"))
    gt_log_set_fp(stderr);
  else {
    *logfp = gt_fa_xfopen(debugfp, "w");
    gt_log_set_fp(*logfp);
  }
}

int gtr_run(GtR *gtr, int argc, const char **argv, GtError *err)
{
  GtToolfunc toolfunc;
  GtTool *tool = NULL;
  char **nargv = NULL;
  void *mem, *map;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(gtr);
  if (gtr->debug)
    enable_logging(gt_str_get(gtr->debugfp), &gtr->logfp);
  gtr->seed = gt_ya_rand_init(gtr->seed);
  gt_log_log("seed=%u", gtr->seed);
  if (gtr->check64bit)
    return check64bit();
  if (gtr->test)
    return run_tests(gtr, err);
  if (gt_str_length(gtr->testspacepeak)) {
    mem = gt_malloc(1 << 26); /* alloc 64 MB */;
    map = gt_fa_xmmap_read(gt_str_get(gtr->testspacepeak), NULL);
    gt_fa_xmunmap(map);
    gt_free(mem);
  }
  if (argc == 0 && !gtr->interactive) {
    gt_error_set(err, "neither tool nor script specified; option -help lists "
                      "possible tools");
    had_err = -1;
  }
  if (!had_err && argc) {
    if (!gtr->tools || !gt_toolbox_has_tool(gtr->tools, argv[0])) {
      /* no tool found -> try to open script */
      if (gt_file_exists(argv[0])) {
        /* run script */
        nargv = gt_cstr_array_prefix_first(argv, gt_error_get_progname(err));
        gt_lua_set_arg(gtr->L, nargv[0], (const char**) nargv+1);
        if (luaL_dofile(gtr->L, argv[0])) {
          /* error */
          gt_assert(lua_isstring(gtr->L, -1)); /* error message on top */
          gt_error_set(err, "could not execute script %s",
                       lua_tostring(gtr->L, -1));
          had_err = -1;
          lua_pop(gtr->L, 1); /* pop error message */
        }
      }
      else {
        /* neither tool nor script found */
        gt_error_set(err, "neither tool nor script '%s' found; option -help "
                          "lists possible tools", argv[0]);
        had_err = -1;
      }
    }
    else {
      /* run tool */
      if (!(toolfunc = gt_toolbox_get(gtr->tools, argv[0]))) {
        tool = gt_toolbox_get_tool(gtr->tools, argv[0]);
        gt_assert(tool);
      }
      nargv = gt_cstr_array_prefix_first(argv, gt_error_get_progname(err));
      gt_error_set_progname(err, nargv[0]);
      if (toolfunc)
        had_err = toolfunc(argc, (const char**) nargv, err);
      else
        had_err = gt_tool_run(tool, argc, (const char**) nargv, err);
    }
  }
  gt_cstr_array_delete(nargv);
  if (!had_err && gtr->interactive) {
    gt_showshortversion(gt_error_get_progname(err));
    gt_lua_set_arg(gtr->L, gt_error_get_progname(err), argv);
    run_interactive_lua_interpreter(gtr->L);
  }
  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

void gtr_delete(GtR *gtr)
{
  if (!gtr) return;
  gt_fa_fclose(gtr->logfp);
  gt_str_delete(gtr->testspacepeak);
  gt_str_delete(gtr->debugfp);
  gt_str_delete(gtr->test_only);
  gt_toolbox_delete(gtr->tools);
  gt_hashmap_delete(gtr->unit_tests);
  if (gtr->L) lua_close(gtr->L);
#ifndef WITHOUT_CAIRO
  gt_style_delete_without_state(gtr->style);
#endif
  gt_free(gtr);
}
