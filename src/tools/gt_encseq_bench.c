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
#include "core/encseq_metadata.h"
#include "core/mathsupport.h"
#include "core/showtime.h"
#include "tools/gt_encseq_bench.h"

typedef struct
{
  unsigned long ccext;
  bool sortlenprepare;
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
                                        "extractions",
                               &arguments->ccext, 0UL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("solepr", "prepare data structure for sequences "
                                         "ordered by their length",
                               &arguments->sortlenprepare, false);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_max_args(op, 1U, 1U);
  return op;
}

static void gt_bench_character_extractions(const GtEncseq *encseq,
                                           unsigned long ccext)
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
  for (idx = 0; idx < ccext; idx++)
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

typedef struct
{
  unsigned long minlength, maxlength, *distrib;
} GtSortedlengthinfo;

void gt_sortedlengthinfo_delete(GtSortedlengthinfo *sortedlengthinfo)
{
  if (sortedlengthinfo != NULL)
  {
    gt_free(sortedlengthinfo->distrib);
    gt_free(sortedlengthinfo);
  }
}

GtSortedlengthinfo *gt_sortedlengthinfo_new(const GtEncseq *encseq,
                                            const char *indexname,
                                            GtError *err)
{
  GtEncseqMetadata *emd;
  GtSortedlengthinfo *sortedlengthinfo = NULL;
  int had_err = 0;

  emd = gt_encseq_metadata_new(indexname, err);
  if (emd == NULL)
  {
    had_err = -1;
  } else
  {
    unsigned long seqlen, previousseqlen = 0, seqnum, numofentries;

    sortedlengthinfo = gt_malloc(sizeof *sortedlengthinfo);
    sortedlengthinfo->minlength = gt_encseq_metadata_min_seq_length(emd);
    sortedlengthinfo->maxlength = gt_encseq_metadata_max_seq_length(emd);
    gt_assert(sortedlengthinfo->minlength > 0 &&
              sortedlengthinfo->minlength <= sortedlengthinfo->maxlength);
    numofentries = sortedlengthinfo->maxlength - sortedlengthinfo->minlength+1;
    sortedlengthinfo->distrib = gt_calloc((size_t) numofentries,
                                          sizeof (*sortedlengthinfo->distrib));
    sortedlengthinfo->distrib -= sortedlengthinfo->minlength;
    for (seqnum = 0; seqnum < gt_encseq_num_of_sequences(encseq); seqnum++)
    {
      seqlen = gt_encseq_seqlength(encseq,seqnum);

      if (seqnum > 0 && previousseqlen > seqlen)
      {
        gt_error_set(err,"sequence %lu of length %lu is longer than "
                         "sequence %lu of length %lu",
                         seqnum-1,previousseqlen,seqnum,seqlen);
        had_err = -1;
        break;
      }
      gt_assert(seqlen >= sortedlengthinfo->minlength &&
                seqlen < numofentries + sortedlengthinfo->minlength);
      sortedlengthinfo->distrib[seqlen]++;
      previousseqlen = seqlen;
    }
    for (seqlen = sortedlengthinfo->minlength;
         seqlen <= sortedlengthinfo->maxlength; seqlen++)
    {
      if (sortedlengthinfo->distrib[seqlen] > 0)
      {
        printf("# %lu %lu\n",seqlen,sortedlengthinfo->distrib[seqlen]);
      }
    }

    sortedlengthinfo->distrib += sortedlengthinfo->minlength;
  }
  gt_encseq_metadata_delete(emd);
  if (had_err)
  {
    gt_sortedlengthinfo_delete(sortedlengthinfo);
    sortedlengthinfo = NULL;
  }
  return sortedlengthinfo;
}

static int gt_encseq_bench_runner(GT_UNUSED int argc, const char **argv,
                                  int parsed_args, void *tool_arguments,
                                  GtError *err)
{
  GtEncseqBenchArguments *arguments = tool_arguments;
  GtEncseqLoader *encseq_loader;
  GtEncseq *encseq;
  int had_err = 0;
  const char *indexname;

  gt_error_check(err);
  gt_assert(arguments != NULL);
  encseq_loader = gt_encseq_loader_new();
  indexname = argv[parsed_args];
  encseq = gt_encseq_loader_load(encseq_loader, indexname, err);
  if (encseq == NULL)
  {
    had_err = -1;
  } else
  {
    if (arguments->sortlenprepare)
    {
      GtSortedlengthinfo *sortedlengthinfo
        = gt_sortedlengthinfo_new(encseq,indexname,err);

      if (sortedlengthinfo == NULL)
      {
        had_err = -1;
      } else
      {
        gt_sortedlengthinfo_delete(sortedlengthinfo);
      }
    }
    if (!had_err && arguments->ccext > 0)
    {
      gt_bench_character_extractions(encseq,arguments->ccext);
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
