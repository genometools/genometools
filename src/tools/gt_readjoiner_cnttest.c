/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/encseq.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "match/esa-seqread.h"
#include "match/rdj-contfind-bottomup.h"
#include "match/rdj-filesuf-def.h"
#include "match/rdj-ovlfind-gusfield.h"
#include "match/rdj-pairwise.h"
#include "match/rdj-cntlist.h"
#include "match/rdj-version.h"
#include "tools/gt_readjoiner_cnttest.h"

typedef enum {
  GT_READJOINER_CNTTEST_SHOWLIST,
  GT_READJOINER_CNTTEST_BRUTEFORCE,
  GT_READJOINER_CNTTEST_KMP,
  GT_READJOINER_CNTTEST_ESA,
} GtReadjoinerCnttestTestspec;

typedef struct {
  GtStr *readset, *teststr;
  GtReadjoinerCnttestTestspec test;
  bool singlestrand;
  bool verbose;
} GtReadjoinerCnttestArguments;

static void* gt_readjoiner_cnttest_arguments_new(void)
{
  GtReadjoinerCnttestArguments *arguments = gt_calloc((size_t)1,
      sizeof *arguments);
  arguments->readset = gt_str_new();
  arguments->teststr = gt_str_new();
  return arguments;
}

static void gt_readjoiner_cnttest_arguments_delete(void *tool_arguments)
{
  GtReadjoinerCnttestArguments *arguments = tool_arguments;
  if (arguments == NULL)
    return;
  gt_str_delete(arguments->readset);
  gt_str_delete(arguments->teststr);
  gt_free(arguments);
}

static GtOptionParser* gt_readjoiner_cnttest_option_parser_new(
    void *tool_arguments)
{
  GtReadjoinerCnttestArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("-readset [option ...]",
      "Test/development tool for readjoiner containments filtering.");

  /* -readset */
  option = gt_option_new_string("readset", "specify the readset name",
      arguments->readset, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -test */
  option = gt_option_new_string("test",
      "select among the available tests:\n"
      "showlist: show content of an cnt list (input: cntlist)\n"
      "bruteforce: memcmp reads vs all suffixes and prefixes (input: encseq)\n"
      "kmp: variant of Knuth-Morris-Pratt (input: encseq)\n"
      "esa: own esa-based algorithm (input: esa)",
      arguments->teststr, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -singlestrand */
  option = gt_option_new_bool("singlestrand", "do not use reads "
      "reverse complements", &arguments->singlestrand, false);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_version_func(op, gt_readjoiner_show_version);
  gt_option_parser_set_max_args(op, 0);

  return op;
}

static int gt_readjoiner_cnttest_arguments_check(GT_UNUSED int rest_argc,
    void *tool_arguments, GT_UNUSED GtError *err)
{
  GtReadjoinerCnttestArguments *arguments = tool_arguments;
  const char *teststr = gt_str_get(arguments->teststr);
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (strcmp(teststr, "showlist") == 0)
    arguments->test = GT_READJOINER_CNTTEST_SHOWLIST;
  else if (strcmp(teststr, "bruteforce") == 0)
    arguments->test = GT_READJOINER_CNTTEST_BRUTEFORCE;
  else if (strcmp(teststr, "kmp") == 0)
    arguments->test = GT_READJOINER_CNTTEST_KMP;
  else if (strcmp(teststr, "esa") == 0)
    arguments->test = GT_READJOINER_CNTTEST_ESA;
  else
  {
    gt_error_set(err, "illegal argument \"%s\" to option -test", teststr);
    had_err = -1;
  }
  return had_err;
}

static int gt_readjoiner_cnttest_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    void *tool_arguments, GT_UNUSED GtError *err)
{
  GtReadjoinerCnttestArguments *arguments = tool_arguments;
  GtEncseqLoader *el = NULL;
  GtEncseq *reads = NULL;
  GtBitsequence *bits = NULL;
  unsigned long nofreads;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (arguments->test == GT_READJOINER_CNTTEST_SHOWLIST)
  {
    GtStr *fn = NULL;
    fn = gt_str_clone(arguments->readset);
    gt_str_append_cstr(fn, GT_READJOINER_SUFFIX_CNTLIST);
    had_err = gt_cntlist_parse(gt_str_get(fn), true, &bits, &nofreads, err);
    gt_str_delete(fn);
  }
  else if (arguments->test == GT_READJOINER_CNTTEST_BRUTEFORCE ||
      arguments->test == GT_READJOINER_CNTTEST_KMP)
  {
    el = gt_encseq_loader_new();
    gt_encseq_loader_drop_description_support(el);
    gt_encseq_loader_disable_autosupport(el);
    if (!arguments->singlestrand)
      gt_encseq_loader_mirror(el);
    reads = gt_encseq_loader_load(el, gt_str_get(arguments->readset), err);
    if (reads == NULL)
      had_err = -1;
    else
    {
      gt_rdj_pairwise_exact(GT_OVLFIND_CNT, reads, !arguments->singlestrand,
          false, arguments->test == GT_READJOINER_CNTTEST_KMP, 1UL, true,
          NULL, NULL, false, NULL, &bits, &nofreads);
    }
    gt_encseq_delete(reads);
    gt_encseq_loader_delete(el);
  }
  else if (arguments->test == GT_READJOINER_CNTTEST_ESA)
  {
    Sequentialsuffixarrayreader *ssar = NULL;
    unsigned long readlength = 0, firstrevcompl = 0;
    GtLogger *verbose_logger = gt_logger_new(arguments->verbose,
        GT_LOGGER_DEFLT_PREFIX, stderr);
    ssar = gt_newSequentialsuffixarrayreaderfromfile(gt_str_get(
          arguments->readset), SARR_LCPTAB | SARR_SUFTAB | SARR_SSPTAB,
        SEQ_scan, verbose_logger, err);
    if (gt_error_is_set(err))
      had_err = -1;
    else
    {
      nofreads = gt_encseq_num_of_sequences(ssar->encseq);
      if (!arguments->singlestrand)
      {
        nofreads = GT_DIV2(nofreads);
        firstrevcompl = nofreads;
      }
      GT_INITBITTAB(bits, nofreads);
      if (!arguments->singlestrand)
      if (gt_encseq_accesstype_get(ssar->encseq) == GT_ACCESS_TYPE_EQUALLENGTH)
        readlength = gt_encseq_seqlength(ssar->encseq, 0);
      (void)gt_contfind_bottomup(ssar, false, bits, arguments->singlestrand ? 0
          : firstrevcompl, readlength);
    }
    if (ssar != NULL)
      gt_freeSequentialsuffixarrayreader(&ssar);
    gt_logger_delete(verbose_logger);
  }
  else
  {
    gt_assert(false);
  }
  if (!had_err)
    had_err = gt_cntlist_show(bits, nofreads, NULL, false, err);
  gt_free(bits);
  return had_err;
}

GtTool* gt_readjoiner_cnttest(void)
{
  return gt_tool_new(gt_readjoiner_cnttest_arguments_new,
                  gt_readjoiner_cnttest_arguments_delete,
                  gt_readjoiner_cnttest_option_parser_new,
                  gt_readjoiner_cnttest_arguments_check,
                  gt_readjoiner_cnttest_runner);
}
