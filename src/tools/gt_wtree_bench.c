/*
  Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>

#include "core/chardef.h"
#include "core/encseq_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/timer_api.h"
#include "core/unused_api.h"
#include "extended/wtree_encseq.h"
#include "tools/gt_wtree_bench.h"

#define WAVELET_BENCH_SIZE 1000000UL
typedef struct {
  GtStr  *safe;
} GtWaveletBenchArguments;

static void* gt_wtree_bench_arguments_new(void)
{
  GtWaveletBenchArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->safe = gt_str_new();
  return arguments;
}

static void gt_wtree_bench_arguments_delete(void *tool_arguments)
{
  GtWaveletBenchArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->safe);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_wtree_bench_option_parser_new(void *tool_arguments)
{
  GtWaveletBenchArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] encoded_sequence",
                            "Testing and benchmarking for wtree.");

  /* -safe */
  option = gt_option_new_string("safe", "safe files to disk, currently not "
                                "implemented",
                                arguments->safe, NULL);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_wtree_bench_arguments_check(int rest_argc,
                                          void *tool_arguments,
                                          GtError *err)
{
  GtWaveletBenchArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (rest_argc == 0) {
    gt_error_set(err, "no encseq basename given");
    had_err = 1;
  }
  else if (rest_argc > 1) {
    gt_error_set(err, "give only one encseq basename");
    had_err = 1;
  }
  if (!had_err && gt_str_length(arguments->safe))
    printf("%s\n", gt_str_get(arguments->safe));

  return had_err;
}

static int gt_wtree_bench_bench_encseq(GtEncseq *es, GtTimer *t,
                                       GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtUword idx, length, pos;
  gt_error_check(err);
  length = gt_encseq_total_length(es);
  gt_timer_start(t);
  for (idx = 0; idx < WAVELET_BENCH_SIZE; idx++) {
    pos = gt_rand_max(length - 1);
    if (gt_encseq_position_is_separator(es, pos, GT_READMODE_FORWARD)) {
      printf("$");
    }
    else {
      printf("%c", gt_encseq_get_decoded_char(es, pos, GT_READMODE_FORWARD));
    }
  }
  printf("\n");
  gt_timer_show_progress_final(t, stderr);
  gt_timer_stop(t);
  return had_err;
}

static int gt_wtree_bench_bench_wtree(GtWtree *wt,
                                      GtError *err,
                                      GtTimer *timer)
{
  int had_err = 0;
  GtUword idx,
                length = gt_wtree_length(wt),
                syms = gt_wtree_num_of_symbols(wt),
                tmp, pos, *max_ranks;
  char c;
  GtWtreeSymbol symbol;
  gt_error_check(err);
  gt_timer_show_progress(timer, "1M random access", stderr);
  printf("\n");
  for (idx = 0; !had_err && idx < WAVELET_BENCH_SIZE; idx++) {
    symbol = gt_wtree_access(wt, gt_rand_max(length-1));
    c = gt_wtree_encseq_unmap_decoded(wt, symbol);
    switch (c) {
      case (char) SEPARATOR:
        printf("$");
        break;
      case (char) UNDEFCHAR:
        gt_error_set(err, "undefined char in sequence, can't print");
        had_err = 1;
        break;
      default:
        printf("%c",c);
    }
  }
  gt_timer_show_progress(timer, "1M random rank", stderr);
  printf("\n");
  for (idx = 0; !had_err && idx < WAVELET_BENCH_SIZE; idx++) {
    symbol = gt_rand_max(syms-1);
    pos = gt_rand_max(length-1);
    tmp = gt_wtree_rank(wt, pos, symbol);
    c = gt_wtree_encseq_unmap_decoded(wt, symbol);
    if (isprint(c))
      printf("rank of %c at "GT_WU": "GT_WU"\n", c, pos, tmp);
    else
      printf("rank of %d at "GT_WU": "GT_WU"\n", c, pos, tmp);
  }
  printf("\n");
  gt_timer_show_progress(timer, "1M random select", stderr);
  max_ranks = gt_malloc((size_t) syms * sizeof (*max_ranks));
  for (idx = 0; !had_err && idx < syms; idx++) {
    max_ranks[idx] = gt_wtree_rank(wt, length - 1, idx);
  }
  printf("\n");
  for (idx = 0; !had_err && idx < WAVELET_BENCH_SIZE; idx++) {
    do {
    symbol = gt_rand_max(syms-1);
    } while (max_ranks[symbol] == 0);
    do {
    pos = gt_rand_max(max_ranks[symbol]);
    } while (pos == 0);
    tmp = gt_wtree_select(wt, pos, symbol);
    c = gt_wtree_encseq_unmap_decoded(wt, symbol);
    if (isprint(c))
      printf("select "GT_WU"th %c: at "GT_WU"\n", pos, c, tmp);
    else
      printf("select "GT_WU"th %d: at "GT_WU"\n", pos, c, tmp);
  }
  printf("\n");
  gt_free(max_ranks);
  return had_err;
}

static int gt_wtree_bench_runner(GT_UNUSED int argc, const char **argv,
                                 int parsed_args,
                                 GT_UNUSED void *tool_arguments,
                                 GT_UNUSED GtError *err)
{
  GT_UNUSED GtWaveletBenchArguments *arguments = tool_arguments;
  int had_err = 0;
  GtEncseq *encseq;
  GtEncseqLoader *el = gt_encseq_loader_new();
  const char *es_basename = argv[parsed_args];
  GtWtree *wt = NULL;
  GtTimer *timer =
    gt_timer_new_with_progress_description("random access encseq 1M");

  gt_error_check(err);
  gt_assert(arguments);

  encseq = gt_encseq_loader_load(el, es_basename, err);
  had_err = gt_wtree_bench_bench_encseq(encseq, timer, err);
  gt_timer_delete(timer);

  if (!had_err) {
    timer = gt_timer_new_with_progress_description("creating wt");
    gt_timer_start(timer);
    wt = gt_wtree_encseq_new(encseq);
    had_err = gt_wtree_bench_bench_wtree(wt, err, timer);
    gt_timer_show_progress_final(timer, stderr);
  }
  gt_timer_delete(timer);
  gt_encseq_delete(encseq);

  gt_encseq_loader_delete(el);
  gt_wtree_delete(wt);

  return had_err;
}

GtTool* gt_wtree_bench(void)
{
  return gt_tool_new(gt_wtree_bench_arguments_new,
                     gt_wtree_bench_arguments_delete,
                     gt_wtree_bench_option_parser_new,
                     gt_wtree_bench_arguments_check,
                     gt_wtree_bench_runner);
}
