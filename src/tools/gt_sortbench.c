/*
  Copyright (c) 2010-2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/qsort_r_api.h"
#include "core/str.h"
#include "core/timer_api.h"
#include "core/radix_sort.h"
#include "core/unused_api.h"
#include "core/qsort-ulong.h"
#include "tools/gt_sortbench.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif

typedef struct {
  GtStr *impl;
  GtUword num_values,
                maxvalue;
  bool use_aqsort,
       use_permute,
       verbose;
} QSortBenchArguments;

static void *gt_sortbench_arguments_new(void)
{
  QSortBenchArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->impl = gt_str_new();
  return arguments;
}

static void gt_sortbench_arguments_delete(void *tool_arguments)
{
  QSortBenchArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->impl);
  gt_free(arguments);
}

static const char *gt_sort_implementation_names[]
    = {"thomas","system","inlinedptr","inlinedarr","direct",
       "radixinplace","radixlsb",NULL};

static GtOptionParser* gt_sortbench_option_parser_new(void *tool_arguments)
{
  QSortBenchArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...]",
                            "Benchmarks quicksort implementations.");

  option = gt_option_new_choice(
                 "impl", "implementation\nchoose from "
                 "thomas|system|inlinedptr|inlinedarr|direct|\n"
                 "radixinplace|radixlsb",
                  arguments->impl,
                  gt_sort_implementation_names[0],
                  gt_sort_implementation_names);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_uword("size", "number of integers to sort",
                               &arguments->num_values, 1000000UL);
  gt_option_parser_add_option(op, option);

  /* default set to ULONG_MAX-1 to facilitate proper testing of the
     radixsort implementations (make sure that higher order bits are set) */
  option = gt_option_new_uword("maxval", "maximal integer to sort",
                               &arguments->maxvalue, ULONG_MAX-1);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("aqsort", "prepare bad input array using the "
                                        "'aqsort' anti-quicksort algorithm",
                               &arguments->use_aqsort, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("permute", "prepare bad input array by "
                                         "permutation of unique items",
                               &arguments->use_permute, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);
  return op;
}

static int gt_sortbench_arguments_check(GT_UNUSED int rest_argc,
                                        void *tool_arguments,
                                        GT_UNUSED GtError *err)
{
  GT_UNUSED QSortBenchArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_assert(arguments != NULL);
  return had_err;
}

static void gt_qbenchsort_permute_ulong_array (GtUword *arr,
                                               GtUword size)
{
  GtUword i, j, c;

  for (i = 0; i < size; ++i)
  {
    j = gt_rand_max(size-1);
    c = arr[i];
    arr[i] = arr[j];
    arr[j] = c;
  }
}

GtUword *val;         /* array, solidified on the fly */
GtUword ncmp;         /* number of comparisons */
GtUword nsolid;       /* number of solid items */
GtUword candidate;    /* pivot candidate */
GtUword gas;          /* gas value = highest sorted value */
#define freeze(x) val[x] = nsolid++

static int gt_sortbench_cmp(const void *px, const void *py)
{
  const GtUword x = *(const GtUword*)px;
  const GtUword y = *(const GtUword*)py;
  ncmp++;
  if (val[x]==gas && val[y]==gas) {
    if (x == candidate) {
      freeze(x);
    } else {
      freeze(y);
    }
  }
  if (val[x] == gas)
    candidate = x;
  else if (val[y] == gas)
    candidate = y;
  if (val[x] < val[y])
  {
    return -1;
  }
  if (val[x] > val[y])
  {
    return 1;
  }
  return 0;
}

static void gt_sortbench_aqsort(GtUword n, GtUword *a)
{
  GtUword i;
  GtUword *ptr = gt_malloc(sizeof (*ptr) * n);

  val = a;
  gas = n-1;
  nsolid = ncmp = candidate = 0;
  for (i=0; i<n; i++) {
    ptr[i] = i;
    val[i] = gas;
  }
  qsort(ptr, (size_t) n, sizeof (*ptr), gt_sortbench_cmp);
  for (i=1UL;i<n;i++)
    gt_assert(val[ptr[i]]==val[ptr[i-1]]+1);
  gt_free(ptr);
}

typedef GtUword Sorttype;

static GtUword cmpcount = 0;

static int qsortcmp(const Sorttype *a,const Sorttype *b,
                    const GT_UNUSED void *data)
{
  cmpcount++;
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

#ifdef  QSORTNAME
#undef  QSORTNAME
#endif

#define QSORTNAME(NAME) sortbench_##NAME

#define sortbench_ARRAY_GET(ARR,RELIDX) ARR[RELIDX]

#define sortbench_ARRAY_SET(ARR,RELIDX,VALUE) ARR[RELIDX] = VALUE

typedef void * QSORTNAME(Datatype);

typedef GtUword QSORTNAME(Sorttype);

static int QSORTNAME(qsortcmparr)(const QSORTNAME(Sorttype) *arr,
                                  GtUword a,
                                  GtUword b,
                                  const GT_UNUSED void *data)
{
  cmpcount++;
  if (arr[a] < arr[b])
  {
    return -1;
  }
  if (arr[a] > arr[b])
  {
    return 1;
  }
  return 0;
}

static void gt_sortbench_verify(GT_UNUSED const GtUword *arr,
                                 GtUword len)
{
  GtUword idx;

  for (idx = 1UL; idx < len; idx++)
  {
    gt_assert(arr[idx-1] <= arr[idx]);
  }
  printf("verified\n");
}

#include "match/qsort-array.gen"

static void check_inlinedarr_qsort(GtUword *arr, GtUword len)
{
  QSORTNAME(gt_inlinedarr_qsort_r) (6UL, false, arr, len, NULL);
  gt_sortbench_verify(arr,len);
}

static void check_direct_qsort(GtUword *arr, GtUword len)
{
  gt_direct_qsort_ulong (6UL, false, arr, len);
  gt_sortbench_verify(arr,len);
}

static int sortcmpwithdata(const void *a,const void *b, GT_UNUSED void *data)
{
  cmpcount++;
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

static void check_thomas_qsort(GtUword *arr, GtUword len)
{
  gt_qsort_r(arr,(size_t) len,sizeof (Sorttype),NULL, sortcmpwithdata);
  gt_sortbench_verify(arr,len);
}

static int sortcmpnodata(const void *a,const void *b)
{
  cmpcount++;
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

static void check_gnu_qsort(GtUword *arr, GtUword len)
{
  qsort(arr,(size_t) len, sizeof (Sorttype), sortcmpnodata);
  gt_sortbench_verify(arr,len);
}

static void check_inlinedptr_qsort(GtUword *arr, GtUword len)
{
  gt_inlined_qsort_r(arr, len, NULL);
  gt_sortbench_verify(arr,len);
}

static void check_radixsort_lsb(GtUword *arr, GtUword len)
{
  /* note that lsb_linear allocates extra temp space of size equal to
     the area to be sorted */
  gt_radixsort_lsb_linear(arr,len);
  gt_sortbench_verify(arr,len);
}

static void check_radixsort_inplace(GtUword *arr, GtUword len)
{
  gt_radixsort_inplace_ulong(arr,len);
  gt_sortbench_verify(arr,len);
}

typedef void (*GtQsortimplementationfunc)(GtUword *,GtUword);

static GtQsortimplementationfunc gt_sort_implementation_funcs[] =
{
  check_thomas_qsort,
  check_gnu_qsort,
  check_inlinedptr_qsort,
  check_inlinedarr_qsort,
  check_direct_qsort,
  check_radixsort_inplace,
  check_radixsort_lsb
};

#define GT_NUM_OF_SORT_IMPLEMENTATIONS\
        (sizeof (gt_sort_implementation_funcs)/\
         sizeof (gt_sort_implementation_funcs[0]))

static int gt_sortbench_runner(GT_UNUSED int argc, GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments, GT_UNUSED GtError *err)
{
  QSortBenchArguments *arguments = tool_arguments;
  int had_err = 0;
  size_t method;
  GtTimer *timer;
  GtUword *array, idx;

  gt_error_check(err);
  gt_assert(arguments);
  if (arguments->verbose) {
    printf("# number of items = "GT_WU"\n", arguments->num_values);
    printf("# max value items = "GT_WU"\n", arguments->maxvalue);
    printf("# implementation = %s\n", gt_str_get(arguments->impl));
  }

  /* prepare benchmark input data */
  array = gt_malloc(sizeof (*array) * arguments->num_values );
  if (arguments->use_aqsort) {
    if (arguments->verbose)
      printf("# using aqsort\n");
    gt_sortbench_aqsort(arguments->num_values, array);
  } else if (arguments->use_permute) {
    if (arguments->verbose)
      printf("# using permuted array\n");
    for (idx = 0; idx < arguments->num_values; idx++) {
      array[idx] = idx;
    }
    gt_qbenchsort_permute_ulong_array(array, arguments->num_values);
  } else {
    if (arguments->verbose)
      printf("# using simple array of random numbers\n");
    /* use seed initialization to make array deterministic */
#ifndef _WIN32
    srand48(366292341);
    for (idx = 0; idx < arguments->num_values; idx++) {
      array[idx] = drand48() * arguments->maxvalue;
    }
#else
    /* XXX: use random instead of drand48() above */
    fprintf(stderr, "drand48() not implemented\n");
    exit(EXIT_FAILURE);
#endif
  }

  timer = gt_timer_new();
  gt_timer_start(timer);
  /* run implementation */
  for (method = 0; method < GT_NUM_OF_SORT_IMPLEMENTATIONS; method++)
  {
    if (strcmp(gt_str_get(arguments->impl),
               gt_sort_implementation_names[method]) == 0)
    {
      gt_sort_implementation_funcs[method](array, arguments->num_values);
      break;
    }
  }
  gt_timer_show(timer, stdout);
  gt_timer_delete(timer);
  gt_free(array);
  if (cmpcount > 0)
  {
    printf("cmpcount = "GT_WU"\n",cmpcount);
  }
  return had_err;
}

GtTool* gt_sortbench(void)
{
  return gt_tool_new(gt_sortbench_arguments_new,
                     gt_sortbench_arguments_delete,
                     gt_sortbench_option_parser_new,
                     gt_sortbench_arguments_check,
                     gt_sortbench_runner);
}
