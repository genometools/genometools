/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"
#include "gtr.h"
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
#include "tools/gt_stat.h"

struct GTR {
  bool test;
  Hashtable *tools,
            *unit_tests;
};

GTR* gtr_new(void)
{
  return xcalloc(1, sizeof(GTR));
}

static void show_tool(void *key, void *value, void *data)
{
  const char *toolname;
  Array *toolnames;
  assert(key && value && data);
  toolname = (const char*) key;
  toolnames = (Array*) data;
  array_add(toolnames, toolname);
}

static void show_option_comments(const char *progname, void *data)
{
  Array *toolnames;
  unsigned long i;
  GTR *gtr;
  assert(data);
  gtr = (GTR*) data;
  toolnames = array_new(sizeof(const char*));
  if (gtr->tools) {
    hashtable_foreach(gtr->tools, show_tool, toolnames);
    printf("\nTools:\n\n");
    assert(array_size(toolnames));
    qsort(array_get_space(toolnames), array_size(toolnames),
          array_elem_size(toolnames), compare);
    for (i = 0; i < array_size(toolnames); i++)
      xputs(*(const char**) array_get(toolnames, i));
  }
  array_free(toolnames);
}

int gtr_parse(GTR *gtr, int argc, char **argv)
{
  int parsed_args;
  OptionParser *op;
  Option *o;
  assert(gtr);
  op = option_parser_new("[option ...] [tool ...] [argument ...]",
                         "The GenomeTools (gt) genome analysis system "
                          "(http://genometools.org).");
  option_parser_set_comment_func(op, show_option_comments, gtr);
  o = option_new_bool("test", "perform unit tests and exit", &gtr->test, false);
  option_parser_add_option(op, o);
  parsed_args = option_parser_parse(op, argc, argv, versionfunc);
  option_parser_free(op);
  return parsed_args - 1;
}

void gtr_register_components(GTR *gtr)
{
  assert(gtr);
  /* add tools */
  hashtable_free(gtr->tools);
  gtr->tools = hashtable_new(HASH_STRING, NULL, NULL);
  hashtable_add(gtr->tools, "bioseq", gt_bioseq);
  hashtable_add(gtr->tools, "cds", gt_cds);
  hashtable_add(gtr->tools, "clean", gt_clean);
  hashtable_add(gtr->tools, "csa", gt_csa);
  hashtable_add(gtr->tools, "eval", gt_eval);
  hashtable_add(gtr->tools, "exercise", gt_exercise);
  hashtable_add(gtr->tools, "extractfeat", gt_extractfeat);
  hashtable_add(gtr->tools, "filter", gt_filter);
  hashtable_add(gtr->tools, "gff3", gt_gff3);
  hashtable_add(gtr->tools, "gtf2gff3", gt_gff3);
  hashtable_add(gtr->tools, "merge", gt_merge);
  hashtable_add(gtr->tools, "mmapandread", gt_mmapandread);
  hashtable_add(gtr->tools, "stat", gt_stat);
  /* add unit tests */
  hashtable_free(gtr->unit_tests);
  gtr->unit_tests = hashtable_new(HASH_STRING, NULL, NULL);
  hashtable_add(gtr->unit_tests, "alignment class", alignment_unit_test);
  hashtable_add(gtr->unit_tests, "array class", array_unit_test);
  hashtable_add(gtr->unit_tests, "bittab class", bittab_unit_test);
  hashtable_add(gtr->unit_tests, "bsearch module", bsearch_unit_test);
  hashtable_add(gtr->unit_tests, "countingsort module", countingsort_unit_test);
  hashtable_add(gtr->unit_tests, "dlist class", dlist_unit_test);
  hashtable_add(gtr->unit_tests, "evaluator class", evaluator_unit_test);
  hashtable_add(gtr->unit_tests, "grep module", grep_unit_test);
  hashtable_add(gtr->unit_tests, "hashtable class", hashtable_unit_test);
  hashtable_add(gtr->unit_tests, "hmm class", hmm_unit_test);
  hashtable_add(gtr->unit_tests, "range class", range_unit_test);
  hashtable_add(gtr->unit_tests, "splicedseq class", splicedseq_unit_test);
  hashtable_add(gtr->unit_tests, "splitter class", splitter_unit_test);
  hashtable_add(gtr->unit_tests, "string class", str_unit_test);
  hashtable_add(gtr->unit_tests, "tokenizer class", tokenizer_unit_test);
}

void run_test(void *key, void *value, void *data)
{
  const char *testname;
  int (*test)(void);
  int rval, *rvalp;
  assert(key && value && data);
  testname = (const char*) key;
  test = value;
  rvalp = (int*) data;
  printf("%s...", testname);
  xfflush(stdout);
  rval = test();
  switch (rval) {
    case EXIT_SUCCESS:
      xputs("ok");
      break;
    case EXIT_FAILURE:
      xputs("error");
      *rvalp = rval;
      break;
    default: assert(0);
  }
  xfflush(stdout);
}

static int run_tests(GTR *gtr)
{
  int rval = EXIT_SUCCESS;
  assert(gtr);

  /* The following type assumptions are made in the GenomeTools library. */
  ensure(sizeof(char) == 1);
  ensure(sizeof(unsigned char) == 1);
  ensure(sizeof(short) == 2);
  ensure(sizeof(unsigned short) == 2);
  ensure(sizeof(int) == 4);
  ensure(sizeof(unsigned int) == 4);
  ensure(sizeof(long) == 4 || sizeof(long) == 8);
  ensure(sizeof(unsigned long) == 4 || sizeof(unsigned long) == 8);
  ensure(sizeof(unsigned long) >= sizeof(size_t));
  ensure(sizeof(long long) == 8);
  ensure(sizeof(unsigned long long) == 8);

  if (gtr->unit_tests) {
    hashtable_foreach(gtr->unit_tests, run_test, &rval);
  }
  return rval;
}

int gtr_run(GTR *gtr, int argc, char **argv)
{
  int (*tool)(int, char**) = NULL;
  char **nargv;
  int rval;
  assert(gtr);
  if (gtr->test) {
    return run_tests(gtr);
  }
  assert(argc);
  if (argc == 1)
    error("no tool specified; option -help lists possible tools");
  if (!gtr->tools || !(tool = hashtable_get(gtr->tools, argv[1])))
    error("tool '%s' not found; option -help lists possible tools", argv[1]);
  assert(argc);
  nargv = cstr_array_prefix_first(argv+1, argv[0]);
  rval = tool(argc-1, nargv);
  cstr_array_free(nargv);
  return rval;
}

void gtr_free(GTR *gtr)
{
  if (!gtr) return;
  hashtable_free(gtr->tools);
  hashtable_free(gtr->unit_tests);
  free(gtr);
}
