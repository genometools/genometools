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
#include "match/rdj-filesuf-def.h"
#include "match/rdj-ovlfind-gusfield.h"
#include "match/rdj-pairwise.h"
#include "match/rdj-spmlist.h"
#include "match/rdj-version.h"
#include "tools/gt_readjoiner_spmtest.h"

typedef enum {
  GT_READJOINER_SPMTEST_SHOWLIST,
  GT_READJOINER_SPMTEST_BRUTEFORCE,
  GT_READJOINER_SPMTEST_KMP,
  GT_READJOINER_SPMTEST_GUSFIELD,
} GtReadjoinerSpmtestTestspec;

typedef struct {
  GtStr *readset, *teststr;
  GtReadjoinerSpmtestTestspec test;
  unsigned int minmatchlength;
  bool singlestrand;
  bool verbose;
} GtReadjoinerSpmtestArguments;

static void* gt_readjoiner_spmtest_arguments_new(void)
{
  GtReadjoinerSpmtestArguments *arguments = gt_calloc((size_t)1,
      sizeof *arguments);
  arguments->readset = gt_str_new();
  arguments->teststr = gt_str_new();
  return arguments;
}

static void gt_readjoiner_spmtest_arguments_delete(void *tool_arguments)
{
  GtReadjoinerSpmtestArguments *arguments = tool_arguments;
  if (arguments == NULL)
    return;
  gt_str_delete(arguments->readset);
  gt_str_delete(arguments->teststr);
  gt_free(arguments);
}

static GtOptionParser* gt_readjoiner_spmtest_option_parser_new(
    void *tool_arguments)
{
  GtReadjoinerSpmtestArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("-readset [option ...]",
      "Test/development tool for readjoiner overlap phase.");

  /* -readset */
  option = gt_option_new_string("readset", "specify the readset name",
      arguments->readset, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -l */
  option = gt_option_new_uint_min("l", "specify the minimum spm length",
      &arguments->minmatchlength, 0, 1U);
  gt_option_parser_add_option(op, option);

  /* -test */
  option = gt_option_new_string("test",
      "select among the available tests:\n"
      "showlist: show content of an spm list (input: spmlist)\n"
      "bruteforce: memcmp all suffixes and prefixes (input: encseq)\n"
      "kmp: variant of Knuth-Morris-Pratt (input: encseq)\n"
      "gusfield: Gusfield all pairs suffix-prefix (input: esa)",
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

static int gt_readjoiner_spmtest_arguments_check(GT_UNUSED int rest_argc,
    void *tool_arguments, GT_UNUSED GtError *err)
{
  GtReadjoinerSpmtestArguments *arguments = tool_arguments;
  const char *teststr = gt_str_get(arguments->teststr);
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (strcmp(teststr, "showlist") == 0)
    arguments->test = GT_READJOINER_SPMTEST_SHOWLIST;
  else if (strcmp(teststr, "bruteforce") == 0)
    arguments->test = GT_READJOINER_SPMTEST_BRUTEFORCE;
  else if (strcmp(teststr, "kmp") == 0)
    arguments->test = GT_READJOINER_SPMTEST_KMP;
  else if (strcmp(teststr, "gusfield") == 0)
    arguments->test = GT_READJOINER_SPMTEST_GUSFIELD;
  else
  {
    gt_error_set(err, "illegal argument \"%s\" to option -test", teststr);
    had_err = -1;
  }
  return had_err;
}

static int gt_readjoiner_spmtest_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    void *tool_arguments, GT_UNUSED GtError *err)
{
  GtReadjoinerSpmtestArguments *arguments = tool_arguments;
  GtEncseqLoader *el = NULL;
  GtEncseq *reads = NULL;
  GtLogger *verbose_logger;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (arguments->test == GT_READJOINER_SPMTEST_SHOWLIST)
  {
    GtStr *fn = NULL;
    fn = gt_str_clone(arguments->readset);
    gt_str_append_cstr(fn, GT_READJOINER_SUFFIX_SPMLIST);
    had_err = gt_spmlist_parse(gt_str_get(fn), 0, gt_spmproc_show_ascii,
        NULL, err);
    gt_str_delete(fn);
    return had_err;
  }
  verbose_logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX,
      stderr);
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
    Sequentialsuffixarrayreader *ssar = NULL;
    switch (arguments->test)
    {
      case GT_READJOINER_SPMTEST_BRUTEFORCE:
        gt_rdj_pairwise_exact(GT_OVLFIND_SPM, reads, !arguments->singlestrand,
            false, false, (unsigned long)arguments->minmatchlength, true,
            gt_spmproc_show_ascii, NULL, false, NULL, NULL, NULL);
        break;
      case GT_READJOINER_SPMTEST_KMP:
        gt_rdj_pairwise_exact(GT_OVLFIND_SPM, reads, !arguments->singlestrand,
            false, true, (unsigned long)arguments->minmatchlength, true,
            gt_spmproc_show_ascii, NULL, false, NULL, NULL, NULL);
        break;
      case GT_READJOINER_SPMTEST_GUSFIELD:
        ssar = gt_newSequentialsuffixarrayreaderfromfile(gt_str_get(
              arguments->readset), SARR_LCPTAB | SARR_SUFTAB | SARR_SSPTAB |
            SARR_ESQTAB, SEQ_mappedboth, verbose_logger, err);
        if (gt_error_is_set(err))
          had_err = -1;
        else
          gt_rdj_gusfield(ssar, (unsigned long)arguments->minmatchlength, true,
              true, arguments->singlestrand ? 0 :
              GT_DIV2(gt_encseq_num_of_sequences(reads)), gt_spmproc_show_ascii,
              NULL);
        if (ssar != NULL)
          gt_freeSequentialsuffixarrayreader(&ssar);
        break;
      default:
        gt_assert(false);
    }
    gt_encseq_delete(reads);
  }
  gt_encseq_loader_delete(el);
  gt_logger_delete(verbose_logger);
  return had_err;
}

GtTool* gt_readjoiner_spmtest(void)
{
  return gt_tool_new(gt_readjoiner_spmtest_arguments_new,
                  gt_readjoiner_spmtest_arguments_delete,
                  gt_readjoiner_spmtest_option_parser_new,
                  gt_readjoiner_spmtest_arguments_check,
                  gt_readjoiner_spmtest_runner);
}
