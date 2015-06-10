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
#include "core/dual-pivot-qsort.h"
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
          runs,
          maxvalue;
  bool use_aqsort,
       use_permute,
       verify,
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
    = {"thomas","system","inlinedptr","inlinedarr","direct","dual-pivot",
       "radixinplace","radixlsb","radixkeypair",NULL};

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
                 "dual-pivot|radixinplace|radixlsb|radixkeypair",
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

  option = gt_option_new_uword("runs",
                           "run sort multiple times as specified by arguments",
                           &arguments->runs, 1UL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("aqsort", "prepare bad input array using the "
                                        "'aqsort' anti-quicksort algorithm",
                               &arguments->use_aqsort, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("permute", "prepare bad input array by "
                                         "permutation of unique items",
                               &arguments->use_permute, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("verify", "verify result by checking order of "
                                        "sorted array",
                               &arguments->verify, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);
  return op;
}

static int gt_sortbench_arguments_check(GT_UNUSED int rest_argc,
                                        void *tool_arguments,
                                        GtError *err)
{
  QSortBenchArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_assert(arguments != NULL);
  if ((arguments->use_aqsort || arguments->use_permute) &&
      strcmp(gt_str_get(arguments->impl),"radixkeypair") == 0)
  {
    gt_error_set(err,"options -aqsort and -permute is not compatible with "
                     "option -impl radixkeypair");
    had_err = -1;
  }
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

/* The following implements the method
   as described in
   @article{MCI:1999,
      author = {McIlroy, M. D.},
      title = {A Killer Adversary for Quicksort},
      journal = {Softw. Pract. Exper.},
      volume = {29},
      number = {4},
      year = {1999},
      issn = {0038-0644},
      pages = {341--344}
     }

   See http://www.cs.dartmouth.edu/~doug/mdmspe.pdf for a preprint. */

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
                                GT_UNUSED const GtUword *arr_copy,
                                GtUword len)
{
  GtUword idx;

  for (idx = 1UL; idx < len; idx++)
  {
    gt_assert(arr[idx-1] <= arr[idx]);
  }
  for (idx = 0UL; idx < len; idx++)
  {
    gt_assert(arr[idx] <= arr_copy[idx]);
  }
  printf("verified\n");
}

#ifndef NDEBUG
typedef uint32_t GtSeedExtendPosition;
typedef uint32_t GtSeedExtendSeqnum;

/* do not change the order of the components */

typedef struct
{
  GtSeedExtendSeqnum bseqnum; /*  2nd important sort criterion */
  GtSeedExtendSeqnum aseqnum; /* most important sort criterion */
  GtSeedExtendPosition bpos;
  GtSeedExtendPosition apos;  /*  3rd important sort criterion */
} GtSeedExtendSeedPair;

/* The following implements the comparison function we want when ordering
   values of type GtSeedExtendSeedPair. */

static int compareseedextendseedpair(const GtSeedExtendSeedPair *p1,
                                     const GtSeedExtendSeedPair *p2)
{
  if (p1->aseqnum < p2->aseqnum)
  {
    return -1;
  }
  if (p1->aseqnum > p2->aseqnum)
  {
    return 1;
  }
  if (p1->bseqnum < p2->bseqnum)
  {
    return -1;
  }
  if (p1->bseqnum > p2->bseqnum)
  {
    return 1;
  }
  if (p1->apos < p2->apos)
  {
    return -1;
  }
  if (p1->apos > p2->apos)
  {
    return 1;
  }
  if (p1->bpos < p2->bpos)
  {
    return -1;
  }
  if (p1->bpos > p2->bpos)
  {
    return 1;
  }
  return 0;
}
#endif

static void gt_sortbench_verify_keypair(GT_UNUSED const Gtuint64keyPair *arr,
                                        GT_UNUSED const Gtuint64keyPair
                                             *arr_copy,
                                        GtUword len)
{
  GtUword idx;

  gt_assert(sizeof (GtSeedExtendSeedPair) == sizeof (*arr));
  for (idx = 1UL; idx < len; idx++)
  {
    /* We check if the order in GtSeedExtendSeedPair is correctly implemented
       by casting a Gtuint64keyPair-pointer to a GtSeedExtendSeedPair-pointer
       and then calling the comparison function. The order on
       Gtuint64keyPair implies the order on GtSeedExtendSeedPair. */
#ifndef NDEBUG
    GtSeedExtendSeedPair *ptr1, *ptr2;
    ptr1 = (GtSeedExtendSeedPair *) (arr + idx - 1);
    ptr2 = (GtSeedExtendSeedPair *) (arr + idx);
#endif
    gt_assert(arr[idx-1].uint64_a < arr[idx].uint64_a ||
              (arr[idx-1].uint64_a == arr[idx].uint64_a &&
               arr[idx-1].uint64_b <= arr[idx].uint64_b));
    gt_assert(compareseedextendseedpair(ptr1,ptr2) <= 0);
  }
  for (idx = 0UL; idx < len; idx++)
  {
    gt_assert(arr[idx].uint64_a == arr_copy[idx].uint64_a);
    gt_assert(arr[idx].uint64_b == arr_copy[idx].uint64_b);
    /*seedptr = (GtSeedExtendSeedPair *) (arr + idx);
    printf("%12u %12u %12u %12u\n",seedptr->bseqnum,seedptr->aseqnum,
                           seedptr->bpos,seedptr->apos);*/
  }
  printf("verified\n");
}

#include "match/qsort-array.gen"

static void run_inlinedarr_qsort(GtUword *arr, GtUword len)
{
  QSORTNAME(gt_inlinedarr_qsort_r) (6UL, false, arr, len, NULL);
}

static void run_direct_qsort(GtUword *arr, GtUword len)
{
  gt_direct_qsort_ulong (6UL, false, arr, len);
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

static void run_thomas_qsort(GtUword *arr, GtUword len)
{
  gt_qsort_r(arr,(size_t) len,sizeof (Sorttype),NULL, sortcmpwithdata);
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

static void run_gnu_qsort(GtUword *arr, GtUword len)
{
  qsort(arr,(size_t) len, sizeof (Sorttype), sortcmpnodata);
}

static void run_inlinedptr_qsort(GtUword *arr, GtUword len)
{
  gt_inlined_qsort_r(arr, len, NULL);
}

typedef void (*GtQsortimplementationfunc)(GtUword *,GtUword);

static GtQsortimplementationfunc gt_sort_implementation_funcs[] =
{
  run_thomas_qsort,
  run_gnu_qsort,
  run_inlinedptr_qsort,
  run_inlinedarr_qsort,
  run_direct_qsort,
  gt_dual_pivot_qsort,
  gt_radixsort_inplace_ulong,
  /* note that gt_radixsort_lsb_linear allocates extra temp space of size
     equal to the area to be sorted */
  gt_radixsort_lsb_linear
};

static int voidqsortcmp(const void *va,const void *vb)
{
  const GtUword a = *((GtUword *) va);
  const GtUword b = *((GtUword *) vb);

  if (a < b)
  {
    return -1;
  }
  if (a > b)
  {
    return 1;
  }
  return 0;
}

static int voidkeypairqsortcmp(const void *va,const void *vb)
{
  const Gtuint64keyPair *a = (const Gtuint64keyPair *) va;
  const Gtuint64keyPair *b = (const Gtuint64keyPair *) vb;

  if (a->uint64_a < b->uint64_a)
  {
    return -1;
  }
  if (a->uint64_a > b->uint64_a)
  {
    return 1;
  }
  if (a->uint64_b < b->uint64_b)
  {
    return -1;
  }
  if (a->uint64_b > b->uint64_b)
  {
    return 1;
  }
  return 0;
}

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
  GtUword *array = NULL, *array_copy = NULL, idx, r;
  Gtuint64keyPair *arraykeypair = NULL, *arraykeypair_copy = NULL;

  gt_error_check(err);
  gt_assert(arguments);
  if (arguments->verbose) {
    printf("# number of items = "GT_WU"\n", arguments->num_values);
    printf("# max value items = "GT_WU"\n", arguments->maxvalue);
    printf("# implementation = %s\n", gt_str_get(arguments->impl));
  }

  /* prepare benchmark input data */
  if (strcmp(gt_str_get(arguments->impl),"radixkeypair") == 0)
  {
    arraykeypair = gt_malloc(sizeof (*arraykeypair) * arguments->num_values);
  } else
  {
    array = gt_malloc(sizeof (*array) * arguments->num_values);
  }
  if (arguments->use_aqsort)
  {
    if (arguments->verbose)
    {
      printf("# using aqsort\n");
    }
    gt_sortbench_aqsort(arguments->num_values, array);
  } else
  {
    if (arguments->use_permute)
    {
      if (arguments->verbose)
      {
        printf("# using permuted array\n");
      }
      gt_assert(array != NULL);
      for (idx = 0; idx < arguments->num_values; idx++)
      {
        array[idx] = idx;
      }
      gt_qbenchsort_permute_ulong_array(array, arguments->num_values);
    } else
    {
      if (arguments->verbose)
      {
        printf("# using simple array of random numbers\n");
      }
      /* use seed initialization to make array deterministic */
#ifndef _WIN32
      srand48(366292341);
      if (arraykeypair != NULL)
      {
        for (idx = 0; idx < arguments->num_values; idx++)
        {
          arraykeypair[idx].uint64_a = drand48() * arguments->maxvalue;
          arraykeypair[idx].uint64_b = drand48() * arguments->maxvalue;
        }
      } else
      {
        gt_assert(array != NULL);
        for (idx = 0; idx < arguments->num_values; idx++)
        {
          array[idx] = drand48() * arguments->maxvalue;
        }
      }
#else
      /* XXX: use random instead of drand48() above */
      fprintf(stderr, "drand48() not implemented\n");
      exit(EXIT_FAILURE);
#endif
    }
  }
  if (arguments->verify)
  {
    if (arraykeypair != NULL)
    {
      arraykeypair_copy = gt_malloc(sizeof (*arraykeypair_copy) *
                                    arguments->num_values);
      for (idx = 0; idx < arguments->num_values; idx++)
      {
        arraykeypair_copy[idx] = arraykeypair[idx];
      }
      qsort(arraykeypair_copy,(size_t) arguments->num_values,
            sizeof *arraykeypair_copy,
            voidkeypairqsortcmp);
    } else
    {
      array_copy = gt_malloc(sizeof (*array_copy) * arguments->num_values);
      gt_assert(array != NULL);
      for (idx = 0; idx < arguments->num_values; idx++)
      {
        array_copy[idx] = array[idx];
      }
      qsort(array_copy,(size_t) arguments->num_values,sizeof *array_copy,
            voidqsortcmp);
    }
  }
  timer = gt_timer_new();
  gt_timer_start(timer);
  /* run implementation */
  if (arraykeypair != NULL)
  {
    for (r = 0; r < arguments->runs; r++)
    {
      gt_radixsort_inplace_Gtuint64keyPair(arraykeypair,arguments->num_values);
    }
  } else
  {
    for (method = 0; method < GT_NUM_OF_SORT_IMPLEMENTATIONS; method++)
    {
      if (strcmp(gt_str_get(arguments->impl),
                 gt_sort_implementation_names[method]) == 0)
      {
        for (r = 0; r < arguments->runs; r++)
        {
          gt_assert(array != NULL && arraykeypair == NULL);
          gt_assert(method < GT_NUM_OF_SORT_IMPLEMENTATIONS);
          gt_sort_implementation_funcs[method](array, arguments->num_values);
        }
        break;
      }
    }
  }
  printf("# TIME %s-t%u-r" GT_WU "-n" GT_WU " overall ",
          gt_str_get(arguments->impl),
#ifdef GT_THREADS_ENABLED
          gt_jobs,
#else
          1U,
#endif
          arguments->runs,arguments->num_values);
  gt_timer_show_formatted(timer,GT_WD ".%02ld\n",stdout);
  gt_timer_delete(timer);
  if (arguments->verify)
  {
    if (arraykeypair != NULL)
    {
      gt_sortbench_verify_keypair(arraykeypair,arraykeypair_copy,
                                  arguments->num_values);
    } else
    {
      gt_sortbench_verify(array,array_copy,arguments->num_values);
    }
  }
  gt_free(array);
  gt_free(array_copy);
  gt_free(arraykeypair);
  gt_free(arraykeypair_copy);
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
