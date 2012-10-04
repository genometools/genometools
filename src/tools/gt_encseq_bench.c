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
#include "core/logger.h"
#include "tools/gt_encseq_bench.h"

typedef struct
{
  unsigned long ccext;
  bool sortlenprepare, verbose;
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

  option = gt_option_new_verbose(&arguments->verbose);
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
  unsigned long minlength, maxlength, numofdifferentseqlen, *seqlenseppos,
                *seqlencount, *seqlentab;
} GtSortedlengthinfo;

void gt_sortedlengthinfo_delete(GtSortedlengthinfo *sortedlengthinfo)
{
  if (sortedlengthinfo != NULL)
  {
    gt_free(sortedlengthinfo->seqlentab);
    gt_free(sortedlengthinfo->seqlenseppos);
    gt_free(sortedlengthinfo->seqlencount);
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
    unsigned long seqlen, previousseqlen = 0, seqnum,
                  countdifferentlength = 1UL, currentpos = 0;

    sortedlengthinfo = gt_malloc(sizeof *sortedlengthinfo);
    sortedlengthinfo->minlength = gt_encseq_metadata_min_seq_length(emd);
    sortedlengthinfo->maxlength = gt_encseq_metadata_max_seq_length(emd);
    gt_assert(sortedlengthinfo->minlength > 0 &&
              sortedlengthinfo->minlength <= sortedlengthinfo->maxlength);
    for (seqnum = 0; seqnum < gt_encseq_num_of_sequences(encseq); seqnum++)
    {
      seqlen = gt_encseq_seqlength(encseq,seqnum);

      if (seqnum > 0)
      {
        if (previousseqlen > seqlen)
        {
          gt_error_set(err,"sequence %lu of length %lu is longer than "
                           "sequence %lu of length %lu",
                           seqnum-1,previousseqlen,seqnum,seqlen);
          had_err = -1;
          break;
        }
        if (previousseqlen < seqlen)
        {
          countdifferentlength++;
        }
      }
      previousseqlen = seqlen;
    }
    if (!had_err)
    {
      sortedlengthinfo->seqlentab
        = gt_calloc((size_t) countdifferentlength,
                    sizeof (*sortedlengthinfo->seqlentab));
    } else
    {
      sortedlengthinfo->seqlentab = NULL;
    }
    if (!had_err && countdifferentlength >= 2UL)
    {
      sortedlengthinfo->seqlenseppos
        = gt_malloc(sizeof (*sortedlengthinfo->seqlenseppos) *
                    (countdifferentlength-1));
      sortedlengthinfo->seqlencount
        = gt_malloc(sizeof (*sortedlengthinfo->seqlencount) *
                    (countdifferentlength-1));
    } else
    {
      sortedlengthinfo->seqlenseppos = NULL;
      sortedlengthinfo->seqlencount = NULL;
    }
    sortedlengthinfo->numofdifferentseqlen = 0;
    if (!had_err)
    {
      gt_assert(sortedlengthinfo->seqlenseppos != NULL &&
                sortedlengthinfo->seqlencount != NULL &&
                sortedlengthinfo->seqlentab != NULL);
      for (seqnum = 0; seqnum < gt_encseq_num_of_sequences(encseq); seqnum++)
      {
        seqlen = gt_encseq_seqlength(encseq,seqnum);

        if (seqnum > 0)
        {
          if (previousseqlen < seqlen)
          {
            gt_assert(sortedlengthinfo->numofdifferentseqlen <
                      countdifferentlength-1);
            sortedlengthinfo->seqlencount[
              sortedlengthinfo->numofdifferentseqlen] = seqnum;
            sortedlengthinfo->seqlenseppos[
              sortedlengthinfo->numofdifferentseqlen++] = currentpos;
            sortedlengthinfo->seqlentab[
               sortedlengthinfo->numofdifferentseqlen] = seqlen;
            gt_assert(gt_encseq_get_encoded_char(encseq,currentpos,
                                                 GT_READMODE_FORWARD)
                      == (GtUchar) SEPARATOR);
          }
          currentpos += 1UL + seqlen;
        } else
        {
          sortedlengthinfo->seqlentab[0] = seqlen;
          currentpos = seqlen;
        }
        previousseqlen = seqlen;
      }
      sortedlengthinfo->numofdifferentseqlen++;
    }
  }
  gt_encseq_metadata_delete(emd);
  if (!had_err)
  {
    return sortedlengthinfo;
  }
  gt_sortedlengthinfo_delete(sortedlengthinfo);
  return NULL;
}

unsigned long gt_sortedlengthinfo_seqnum(const GtEncseq *encseq,
                                   const GtSortedlengthinfo *sortedlengthinfo,
                                   unsigned long position)
{
  unsigned long recordnum, totallength = gt_encseq_total_length(encseq);

  gt_assert(position < totallength);
  if (gt_encseq_get_encoded_char(encseq,position,GT_READMODE_FORWARD)
                    == (GtUchar) SEPARATOR)
  {
    return ULONG_MAX;
  }
  recordnum = gt_encseq_sep2seqnum(sortedlengthinfo->seqlenseppos,
                                   sortedlengthinfo->numofdifferentseqlen,
                                   totallength,
                                   position);
  if (recordnum == 0)
  {
    return position/(sortedlengthinfo->seqlentab[0]+1);
  } else
  {
    return sortedlengthinfo->seqlencount[recordnum-1] +
           (position - sortedlengthinfo->seqlenseppos[recordnum-1])/
           (sortedlengthinfo->seqlentab[recordnum]+1);
  }
}

unsigned long gt_sortedlengthinfo_seqstart(const GtEncseq *encseq,
                                   const GtSortedlengthinfo *sortedlengthinfo,
                                   unsigned long seqnum)
{
  if (seqnum == 0)
  {
    return 0;
  } else
  {
    unsigned long numberofsequences = gt_encseq_num_of_sequences(encseq),
                  recordnum;

    gt_assert(seqnum < numberofsequences);
    recordnum = gt_encseq_sep2seqnum(sortedlengthinfo->seqlencount,
                                     sortedlengthinfo->numofdifferentseqlen,
                                     numberofsequences,
                                     seqnum);
    if (recordnum == 0)
    {
      return (sortedlengthinfo->seqlentab[0] + 1) * seqnum;
    }
    gt_assert(seqnum >= sortedlengthinfo->seqlencount[recordnum-1]);
    return sortedlengthinfo->seqlenseppos[recordnum-1] +
           (sortedlengthinfo->seqlentab[recordnum] + 1) *
           (seqnum - sortedlengthinfo->seqlencount[recordnum-1]) + 1;
  }
}

static int gt_encseq_bench_runner(GT_UNUSED int argc, const char **argv,
                                  int parsed_args, void *tool_arguments,
                                  GtError *err)
{
  GtEncseqBenchArguments *arguments = tool_arguments;
  GtEncseqLoader *encseq_loader;
  GtEncseq *encseq;
  GtLogger *logger = NULL;
  const char *indexname;
  int had_err = 0;

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
      unsigned long seqnum, position;
      GtSortedlengthinfo *sortedlengthinfo
        = gt_sortedlengthinfo_new(encseq,indexname,err);

      if (sortedlengthinfo == NULL)
      {
        had_err = -1;
      } else
      {
        logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX,
                               stdout);
      }
      gt_logger_log(logger,"prepare sorted length info");
      for (position = 0; !had_err && position < gt_encseq_total_length(encseq);
           position++)
      {
        seqnum = gt_sortedlengthinfo_seqnum(encseq,sortedlengthinfo,position);
        if (seqnum != ULONG_MAX)
        {
          unsigned long seqnum2 = gt_encseq_seqnum(encseq,position);
          if (seqnum2 != seqnum)
          {
            fprintf(stderr,"position %lu: seqnum2 = %lu != %lu = seqnum\n",
                    position,seqnum,seqnum2);
            exit(GT_EXIT_PROGRAMMING_ERROR);
          }
        }
      }
      gt_logger_log(logger,"perform seqnum check");
      for (seqnum = 0; !had_err && seqnum < gt_encseq_num_of_sequences(encseq);
           seqnum++)
      {
        unsigned long seqstart = gt_encseq_seqstartpos(encseq,seqnum);
        unsigned long seqstart2 = gt_sortedlengthinfo_seqstart(encseq,
                                   sortedlengthinfo,seqnum);
        if (seqstart != seqstart2)
        {
          fprintf(stderr,"seqnum=%lu: seqstart = %lu != %lu = seqstart2\n",
                          seqnum,seqstart,seqstart2);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
      }
      gt_logger_log(logger,"perform seqstart check");
      gt_sortedlengthinfo_delete(sortedlengthinfo);
    }
    if (!had_err && arguments->ccext > 0)
    {
      gt_logger_log(logger,"perform character extractions");
      gt_bench_character_extractions(encseq,arguments->ccext);
    }
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(encseq_loader);
  gt_logger_delete(logger);
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
