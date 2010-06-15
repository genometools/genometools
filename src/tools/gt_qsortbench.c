/*
  Copyright (c) CCYY YOUR NAME HERE <user@your.dom.ain>
  Copyright (c) CCYY Center for Bioinformatics, University of Hamburg

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
#include "core/qsort_r.h"
#include "core/ma.h"
#include "core/str.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "tools/gt_qsortbench.h"

typedef struct {
  GtStr *impl;
} QSortBenchArguments;

static void* gt_qsortbench_arguments_new(void)
{
  QSortBenchArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->impl = gt_str_new();
  return arguments;
}

static void gt_qsortbench_arguments_delete(void *tool_arguments)
{
  QSortBenchArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->impl);
  gt_free(arguments);
}

static GtOptionParser* gt_qsortbench_option_parser_new(void *tool_arguments)
{
  QSortBenchArguments *arguments = tool_arguments;
  GtOptionParser *op;
  static const char *inputs[] = {
    "thomas",
    "inlined",
    "system",
    NULL
  };
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...]", /* XXX */
                            "Benchmarks quicksort implementations."); /* XXX */

  option = gt_option_new_choice("impl", "implementation\n"
                                       "choose from thomas|system|inlined",
                             arguments->impl, inputs[0], inputs);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_qsortbench_arguments_check(GT_UNUSED int rest_argc,
                                       GT_UNUSED void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GT_UNUSED QSortBenchArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
}

typedef unsigned long Sorttype;

static int qsortcmp (const Sorttype *a,const Sorttype *b,
                     const GT_UNUSED void *data)
{
  if ((*a) < (*b))
  {
    return -1;
  }
  if ((*a) > (*b))
  {
    return 1;
  }
  return 0;
}

static int qsortcmpnodata (const void *a,const void *b)
{
  if (*(Sorttype*)a < *((Sorttype*)b))
  {
    return -1;
  }
  if (*(Sorttype*)a > *(Sorttype*)b)
  {
    return 1;
  }
  return 0;
}

static int qsortcmpwithdata (const void *a,const void *b, GT_UNUSED void *data)
{
  if (*(Sorttype*)a < *((Sorttype*)b))
  {
    return -1;
  }
  if (*(Sorttype*)a > *(Sorttype*)b)
  {
    return 1;
  }
  return 0;
}

#include "match/qsort-inplace.gen"
#define MAXSIZE 1000000

static void check_inlined_qsort(void)
{
  int times;
  unsigned long idx, n = (unsigned long) MAXSIZE, a[MAXSIZE];

  for (times = 0; times < 100; times++)
  {
    for (idx = 0; idx < n; idx++)
    {
      a[idx] = gt_rand_max(1000UL);
    }
    gt_inlined_qsort_r (a,n,NULL);
    for (idx = 1UL; idx < n; idx++)
    {
      gt_assert(a[idx-1] <= a[idx]);
    }
  }
}

static void check_thomas_qsort(void)
{
  int times;
  unsigned long idx, n = (unsigned long) MAXSIZE, a[MAXSIZE];

  for (times = 0; times < 100; times++)
  {
    for (idx = 0; idx < n; idx++)
    {
      a[idx] = gt_rand_max(1000UL);
    }
    gt_qsort_r(a,n,sizeof (Sorttype),NULL, qsortcmpwithdata);
    for (idx = 1UL; idx < n; idx++)
    {
      gt_assert(a[idx-1] <= a[idx]);
    }
  }
}

static void check_gnu_qsort(void)
{
  int times;
  unsigned long idx, n = (unsigned long) MAXSIZE, a[MAXSIZE];

  for (times = 0; times < 100; times++)
  {
    for (idx = 0; idx < n; idx++)
    {
      a[idx] = gt_rand_max(1000UL);
    }
    qsort (a,n, sizeof (Sorttype), qsortcmpnodata);
    for (idx = 1UL; idx < n; idx++)
    {
      gt_assert(a[idx-1] <= a[idx]);
    }
  }
}

static int gt_qsortbench_runner(GT_UNUSED int argc, GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments, GT_UNUSED GtError *err)
{
  QSortBenchArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (strcmp(gt_str_get(arguments->impl), "thomas") == 0) {
    check_thomas_qsort();
  }
  if (strcmp(gt_str_get(arguments->impl), "system") == 0) {
    check_gnu_qsort();
  }
  if (strcmp(gt_str_get(arguments->impl), "inlined") == 0) {
    check_inlined_qsort();
  }

  return had_err;
}

GtTool* gt_qsortbench(void)
{
  return gt_tool_new(gt_qsortbench_arguments_new,
                  gt_qsortbench_arguments_delete,
                  gt_qsortbench_option_parser_new,
                  gt_qsortbench_arguments_check,
                  gt_qsortbench_runner);
}
