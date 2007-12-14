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
#include "libgtcore/array2dim.h"
#include "libgtcore/bitpackarray.h"
#include "libgtcore/bitpackstring.h"
#include "libgtcore/bittab.h"
#include "libgtcore/bsearch.h"
#include "libgtcore/countingsort.h"
#include "libgtcore/cstr.h"
#include "libgtcore/discdistri.h"
#include "libgtcore/dlist.h"
#include "libgtcore/dynbittab.h"
#include "libgtcore/ensure.h"
#include "libgtcore/fa.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/getbasename.h"
#include "libgtcore/gtdatapath.h"
#include "libgtcore/grep.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/log.h"
#include "libgtcore/range.h"
#include "libgtcore/safearith.h"
#include "libgtcore/splitter.h"
#include "libgtcore/tokenizer.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xansi.h"
#include "libgtext/alignment.h"
#include "libgtext/evaluator.h"
#include "libgtext/genome_node_iterator.h"
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
#include "tools/gt_extractseq.h"
#include "tools/gt_filter.h"
#include "tools/gt_gff3.h"
#include "tools/gt_gff3_to_gtf.h"
#include "tools/gt_gtf_to_gff3.h"
#include "tools/gt_ltrharvest.h"
#include "tools/gt_merge.h"
#include "tools/gt_mkfmindex.h"
#include "tools/gt_mmapandread.h"
#include "tools/gt_mutate.h"
#include "tools/gt_splitfasta.h"
#include "tools/gt_splicesiteinfo.h"
#include "tools/gt_stat.h"
#include "tools/gt_suffixerator.h"
#include "tools/gt_packedindex.h"
#include "tools/gt_uniq.h"
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
#ifdef LIBGTVIEW
  Str *config_file = NULL;
  int had_err = 0;
#endif
  gtr = ma_calloc(1, sizeof (GTR));
  gtr->testspacepeak = str_new();
  gtr->L = luaL_newstate();
  assert(gtr->L); /* XXX: proper error message  */
  luaL_openlibs(gtr->L); /* open the standard libraries */
  luaopen_gt(gtr->L); /* open all GenomeTools libraries */
  luaopen_lfs(gtr->L); /* open Lua filesystem */
#ifdef LIBGTVIEW
  if (!(gtr->config = config_new_with_state(gtr->L, err)))
    had_err = -1;
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
        put_config_in_registry(gtr->L, gtr->config);
    }
  }
  str_delete(config_file);
  if (had_err) {
    ma_free(gtr);
    return NULL;
  }
#endif
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
  gtr->tools = toolbox_new();
  toolbox_add(gtr->tools, "bioseq", gt_bioseq);
  toolbox_add(gtr->tools, "cds", gt_cds);
  toolbox_add(gtr->tools, "chseqids", gt_chseqids);
  toolbox_add(gtr->tools, "clean", gt_clean);
  toolbox_add(gtr->tools, "csa", gt_csa);
  toolbox_add(gtr->tools, "dev", gt_dev);
  toolbox_add(gtr->tools, "eval", gt_eval);
  toolbox_add(gtr->tools, "exercise", gt_exercise);
  toolbox_add(gtr->tools, "extractfeat", gt_extractfeat);
  toolbox_add(gtr->tools, "extractseq", gt_extractseq);
  toolbox_add(gtr->tools, "filter", gt_filter);
  toolbox_add(gtr->tools, "gff3", gt_gff3);
  toolbox_add(gtr->tools, "gff3_to_gtf", gt_gff3_to_gtf);
  toolbox_add(gtr->tools, "gtf_to_gff3", gt_gtf_to_gff3);
  toolbox_add(gtr->tools, "ltrharvest", gt_ltrharvest);
  toolbox_add(gtr->tools, "merge", gt_merge);
  toolbox_add(gtr->tools, "mmapandread", gt_mmapandread);
  toolbox_add(gtr->tools, "mutate", gt_mutate);
  toolbox_add(gtr->tools, "splitfasta", gt_splitfasta);
  toolbox_add(gtr->tools, "splicesiteinfo", gt_splicesiteinfo);
  toolbox_add(gtr->tools, "stat", gt_stat);
  toolbox_add(gtr->tools, "suffixerator", gt_suffixerator);
  toolbox_add(gtr->tools, "packedindex", gt_packedindex);
  toolbox_add(gtr->tools, "mkfmindex", gt_mkfmindex);
  toolbox_add(gtr->tools, "uniq", gt_uniq);
  toolbox_add(gtr->tools, "uniquesub", gt_uniquesub);
#ifdef LIBGTVIEW
  toolbox_add(gtr->tools, "view", gt_view);
#endif
  /* add unit tests */
  hashtable_delete(gtr->unit_tests);
  gtr->unit_tests = hashtable_new(HASH_STRING, NULL, NULL);
  hashtable_add(gtr->unit_tests, "alignment class", alignment_unit_test);
  hashtable_add(gtr->unit_tests, "array class", array_unit_test);
  hashtable_add(gtr->unit_tests, "array example", array_example);
  hashtable_add(gtr->unit_tests, "array2dim example", array2dim_example);
  hashtable_add(gtr->unit_tests, "bit pack array class",
                bitPackArray_unit_test);
  hashtable_add(gtr->unit_tests, "bit pack string module",
                bitPackString_unit_test);
  hashtable_add(gtr->unit_tests, "bittab class", bittab_unit_test);
  hashtable_add(gtr->unit_tests, "bittab example", bittab_example);
  hashtable_add(gtr->unit_tests, "bsearch module", bsearch_unit_test);
  hashtable_add(gtr->unit_tests, "countingsort module", countingsort_unit_test);
  hashtable_add(gtr->unit_tests, "disc distri class", discdistri_unit_test);
  hashtable_add(gtr->unit_tests, "dlist class", dlist_unit_test);
  hashtable_add(gtr->unit_tests, "dlist example", dlist_example);
  hashtable_add(gtr->unit_tests, "dynamic bittab class", dynbittab_unit_test);
  hashtable_add(gtr->unit_tests, "evaluator class", evaluator_unit_test);
  hashtable_add(gtr->unit_tests, "genome node iterator example",
                genome_node_iterator_example);
  hashtable_add(gtr->unit_tests, "getbasename module", getbasename_unit_test);
  hashtable_add(gtr->unit_tests, "grep module", grep_unit_test);
  hashtable_add(gtr->unit_tests, "hashtable class", hashtable_unit_test);
  hashtable_add(gtr->unit_tests, "hmm class", hmm_unit_test);
  hashtable_add(gtr->unit_tests, "range class", range_unit_test);
  hashtable_add(gtr->unit_tests, "safearith module", safearith_unit_test);
  hashtable_add(gtr->unit_tests, "safearith example", safearith_example);
  hashtable_add(gtr->unit_tests, "splicedseq class", splicedseq_unit_test);
  hashtable_add(gtr->unit_tests, "splitter class", splitter_unit_test);
  hashtable_add(gtr->unit_tests, "string class", str_unit_test);
  hashtable_add(gtr->unit_tests, "tokenizer class", tokenizer_unit_test);
#ifdef LIBGTVIEW
  hashtable_add(gtr->unit_tests, "block class", block_unit_test);
  hashtable_add(gtr->unit_tests, "config class", config_unit_test);
  hashtable_add(gtr->unit_tests, "diagram class", diagram_unit_test);
  hashtable_add(gtr->unit_tests, "element class", element_unit_test);
  hashtable_add(gtr->unit_tests, "feature index class",
                feature_index_unit_test);
  hashtable_add(gtr->unit_tests, "line class", line_unit_test);
  hashtable_add(gtr->unit_tests, "track class", track_unit_test);
#endif
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
  Tool tool = NULL;
  char **nargv = NULL;
  void *mem, *map;
  int had_err = 0;
  error_check(err);
  assert(gtr);
  if (gtr->debug)
    log_enable();
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
    if (!gtr->tools || !(tool = toolbox_get(gtr->tools, argv[0]))) {
      /* no tool found -> try to open script */
      if (file_exists(argv[0])) {
        /* run script */
        nargv = cstr_array_prefix_first(argv, error_get_progname(err));
        set_arg_in_lua_interpreter(gtr->L, nargv[0], (const char**) nargv+1);
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
      nargv = cstr_array_prefix_first(argv, error_get_progname(err));
      error_set_progname(err, nargv[0]);
      had_err = tool(argc, (const char**) nargv, err);
    }
  }
  cstr_array_delete(nargv);
  if (!had_err && gtr->interactive) {
    showshortversion(error_get_progname(err));
    set_arg_in_lua_interpreter(gtr->L, error_get_progname(err), argv);
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
