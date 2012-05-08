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
#include "core/fasta_separator.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "match/sfx-bentsedg.h"
#include "match/sfx-suffixgetset.h"
#include "match/sfx-strategy.h"
#include "tools/gt_seqorder.h"

typedef struct {
  bool invert, sort, revsort, shuffle;
} GtSeqorderArguments;

static void* gt_seqorder_arguments_new(void)
{
  GtSeqorderArguments *arguments = gt_calloc((size_t)1, sizeof *arguments);
  return arguments;
}

static void gt_seqorder_arguments_delete(void *tool_arguments)
{
  GtSeqorderArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_free(arguments);
}

static GtOptionParser* gt_seqorder_option_parser_new(void *tool_arguments)
{
  GtSeqorderArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *invert_option, *sort_option, *revsort_option, *shuffle_option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("(-invert|-sort|-revsort|-shuffle) encseq",
                            "Output sequences as MultiFasta in specified "
                            "order.");

  /* -invert */
  invert_option = gt_option_new_bool("invert", "invert order of sequences",
                           &arguments->invert, false);
  gt_option_parser_add_option(op, invert_option);

  /* -sort */
  sort_option = gt_option_new_bool("sort", "sort sequences lexicographically",
                           &arguments->sort, false);
  gt_option_exclude(sort_option, invert_option);
  gt_option_parser_add_option(op, sort_option);

  /* -revsort */
  revsort_option = gt_option_new_bool("revsort", "sort sequences in reverse "
                           "lexicographic order", &arguments->revsort, false);
  gt_option_exclude(revsort_option, invert_option);
  gt_option_exclude(revsort_option, sort_option);
  gt_option_parser_add_option(op, revsort_option);

  /* -shuffle */
  shuffle_option = gt_option_new_bool("shuffle", "shuffle sequences "
                           "pseudo-randomly", &arguments->shuffle, false);
  gt_option_exclude(shuffle_option, invert_option);
  gt_option_exclude(shuffle_option, sort_option);
  gt_option_exclude(shuffle_option, revsort_option);
  gt_option_parser_add_option(op, shuffle_option);

  gt_option_parser_set_min_max_args(op, 1U, 1U);

  return op;
}

static int gt_seqorder_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtSeqorderArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments != NULL);

  if (!(arguments->invert || arguments->sort || arguments->revsort ||
      arguments->shuffle))
  {
    had_err = 1;
    gt_error_set(err, "order option needed: -invert|-sort|-revsort|-shuffle");
  }

  return had_err;
}

static void gt_seqorder_sort(GtSuffixsortspace *suffixsortspace,
    GtEncseq *encseq)
{
  unsigned long i;
  Sfxstrategy sfxstrategy;

  defaultsfxstrategy(&sfxstrategy, false);
  for (i = 0; i < gt_encseq_num_of_sequences(encseq); i++)
    gt_suffixsortspace_setdirect(suffixsortspace, i,
        gt_encseq_seqstartpos(encseq, i));
  gt_sortallsuffixesfromstart(suffixsortspace,
      gt_encseq_num_of_sequences(encseq), encseq, GT_READMODE_FORWARD, NULL, 0,
      &sfxstrategy, NULL, NULL, NULL);
}

static void gt_seqorder_get_shuffled_seqnums(unsigned long nofseqs,
    unsigned long *seqnums)
{
  unsigned long i, j;

  gt_assert(seqnums != NULL);
  seqnums[0] = 0;
  for (i = 1UL; i < nofseqs; i++)
  {
    j = gt_rand_max(i);
    seqnums[i] = seqnums[j];
    seqnums[j] = i;
  }
}

static void gt_seqorder_output(unsigned long seqnum, GtEncseq *encseq)
{
  GtEncseqReader *esr;
  unsigned long startpos, len, desclen = 0;
  const char *desc = NULL;
  unsigned long i;

  startpos = gt_encseq_seqstartpos(encseq, seqnum);
  len = gt_encseq_seqlength(encseq, seqnum);
  gt_xfputc(GT_FASTA_SEPARATOR, stdout);
  if (gt_encseq_has_description_support(encseq))
  {
    desc = gt_encseq_description(encseq, &desclen, seqnum);
    gt_xfwrite(desc, (size_t)1, (size_t)desclen, stdout);
  }
  gt_xfputc('\n', stdout);
  esr = gt_encseq_create_reader_with_readmode(encseq, GT_READMODE_FORWARD,
      startpos);
  for (i = 0; i < len; i++)
  {
    gt_xfputc(gt_encseq_reader_next_decoded_char(esr), stdout);
  }
  gt_encseq_reader_delete(esr);
  gt_xfputc('\n', stdout);
}

static int gt_seqorder_runner(GT_UNUSED int argc, const char **argv,
    int parsed_args, void *tool_arguments, GtError *err)
{
  GtSeqorderArguments *arguments = tool_arguments;
  int had_err = 0;
  GtEncseq *encseq;
  GtEncseqLoader *loader;
  unsigned long i, nofseqs;

  gt_error_check(err);
  gt_assert(arguments != NULL);

  /* load encseq */
  loader = gt_encseq_loader_new();
  encseq = gt_encseq_loader_load(loader, argv[parsed_args], err);
  if (encseq == NULL)
    had_err = -1;
  if (had_err == 0 && !gt_encseq_has_description_support(encseq))
    gt_warning("%s has no description support", argv[parsed_args]);
  if (!had_err)
  {
    nofseqs = gt_encseq_num_of_sequences(encseq);
    if (arguments->invert)
    {
      for (i = nofseqs; i > 0; i--)
        gt_seqorder_output(i - 1, encseq);
    }
    else if (arguments->shuffle)
    {
      unsigned long *seqnums;
      seqnums = gt_malloc(sizeof (unsigned long) * nofseqs);
      gt_seqorder_get_shuffled_seqnums(nofseqs, seqnums);
      for (i = 0; i < nofseqs; i++)
        gt_seqorder_output(seqnums[i], encseq);
      gt_free(seqnums);
    }
    else
    {
      GtSuffixsortspace *suffixsortspace;
      gt_assert(arguments->sort || arguments->revsort);
      suffixsortspace
        = gt_suffixsortspace_new(nofseqs,
                                 /* Use iterator over sequence separators:
                                    saves a lot of binary searches */
                                 gt_encseq_seqstartpos(encseq, nofseqs-1),
                                 false,NULL);
      gt_seqorder_sort(suffixsortspace, encseq);
      if (arguments->sort)
        for (i = 0; i < nofseqs; i++)
          gt_seqorder_output(gt_encseq_seqnum(encseq,
                gt_suffixsortspace_getdirect(suffixsortspace, i)), encseq);
      else
        for (i = nofseqs; i > 0; i--)
          gt_seqorder_output(gt_encseq_seqnum(encseq,
                gt_suffixsortspace_getdirect(suffixsortspace, i - 1)), encseq);
      gt_suffixsortspace_delete(suffixsortspace, false);
    }
  }

  gt_encseq_loader_delete(loader);
  gt_encseq_delete(encseq);
  return had_err;
}

GtTool* gt_seqorder(void)
{
  return gt_tool_new(gt_seqorder_arguments_new,
                  gt_seqorder_arguments_delete,
                  gt_seqorder_option_parser_new,
                  gt_seqorder_arguments_check,
                  gt_seqorder_runner);
}
