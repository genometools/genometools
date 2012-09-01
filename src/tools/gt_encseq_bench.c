/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/unused_api.h"
#include "core/encseq.h"
#include "core/mathsupport.h"
#include "core/showtime.h"
#include "tools/gt_encseq_bench.h"

typedef struct
{
  unsigned long ccext;
} GtEncseqBenchArguments;

static void* gt_encseq_bench_arguments_new(void)
{
  GtEncseqBenchArguments *arguments = gt_malloc(sizeof *arguments);
  return arguments;
}

static void gt_encseq_bench_arguments_delete(void *tool_arguments)
{
  GtEncseqBenchArguments *arguments = tool_arguments;

  if (arguments != NULL)
  {
    gt_free(arguments);
  }
}

static GtOptionParser* gt_encseq_bench_option_parser_new(void *tool_arguments)
{
  GtEncseqBenchArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] indexname",
                         "Perform benchmark on extractions from encseq.");

  option = gt_option_new_ulong("ccext", "specify number of random character "
                                        " extractions",
                               &arguments->ccext, 0UL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  gt_option_parser_set_min_max_args(op, 1U, 1U);
  return op;
}

static int gt_encseq_bench_runner(GT_UNUSED int argc, const char **argv,
                                  int parsed_args, void *tool_arguments,
                                  GT_UNUSED GtError *err)
{
  GtEncseqBenchArguments *arguments = tool_arguments;
  int had_err = 0;
  GtEncseqLoader *encseq_loader;
  GtEncseq *encseq;

  gt_error_check(err);
  gt_assert(arguments != NULL);

  encseq_loader = gt_encseq_loader_new();
  encseq = gt_encseq_loader_load(encseq_loader, argv[parsed_args], err);
  if (encseq == NULL)
  {
    had_err = -1;
  } else
  {
    unsigned long idx, ccsum = 0,
                  totallength = gt_encseq_total_length(encseq);
    GtTimer *timer = NULL;

    if (gt_showtime_enabled())
    {
      timer = gt_timer_new_with_progress_description("run character "
                                                     "extractions");
      gt_timer_start(timer);
    }
    for (idx = 0; idx < arguments->ccext; idx++)
    {
      GtUchar cc;
      unsigned long pos = gt_rand_max(totallength-1);

      cc = gt_encseq_get_encoded_char(encseq,pos,GT_READMODE_FORWARD);
      ccsum += (unsigned long) cc;
    }
    printf("ccsum=%lu\n",ccsum);
    if (timer != NULL)
    {
      gt_timer_show_progress_final(timer, stdout);
      gt_timer_delete(timer);
    }
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(encseq_loader);
  return had_err;
}

GtTool* gt_encseq_bench(void)
{
  return gt_tool_new(gt_encseq_bench_arguments_new,
                  gt_encseq_bench_arguments_delete,
                  gt_encseq_bench_option_parser_new,
                  NULL,
                  gt_encseq_bench_runner);
}
