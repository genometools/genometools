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
#include "core/qsort_r_api.h"
#include "core/ma_api.h"
#include "core/mathsupport_api.h"
#include "core/minmax_api.h"
#include "core/parseutils.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "match/sfx-bentsedg.h"
#include "match/sfx-suffixgetset.h"
#include "match/sfx-strategy.h"
#include "tools/gt_seqorder.h"

typedef struct {
  bool invert, sort, revsort, shuffle, sorthdr, sorthdrnum,
       sort_length;
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
  GtOption *invert_option, *sort_option, *sorthdr_option, *sorthdrnum_option,
           *revsort_option, *shuffle_option, *sort_length_option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("(-invert|-sort|-revsort|-shuffle|-sorthdr"
                            "|-sorthdrnum) encseq",
                            "Output sequences as MultiFasta in specified "
                            "order.");

  /* -invert */
  invert_option = gt_option_new_bool("invert", "invert order of sequences",
                                     &arguments->invert, false);
  gt_option_parser_add_option(op, invert_option);

  /* -sort */
  sort_option = gt_option_new_bool("sort",
                                   "sort sequences lexicographically "
                                   "(by actual sequence)",
                                   &arguments->sort, false);
  gt_option_exclude(sort_option, invert_option);
  gt_option_parser_add_option(op, sort_option);

  /* -revsort */
  revsort_option = gt_option_new_bool("revsort", "sort sequences in reverse "
                                      "lexicographic order",
                                      &arguments->revsort, false);
  gt_option_exclude(revsort_option, invert_option);
  gt_option_exclude(revsort_option, sort_option);
  gt_option_parser_add_option(op, revsort_option);

  /* -sorthdr */
  sorthdr_option = gt_option_new_bool("sorthdr",
                                      "sort sequences lexicographically "
                                      "by sequence header",
                                      &arguments->sorthdr, false);
  gt_option_exclude(sorthdr_option, invert_option);
  gt_option_exclude(sorthdr_option, sort_option);
  gt_option_exclude(sorthdr_option, revsort_option);
  gt_option_parser_add_option(op, sorthdr_option);

    /* -sorthdrnum */
  sorthdrnum_option = gt_option_new_bool("sorthdrnum",
                                         "sort sequences numerically "
                                         "by sequence header",
                                         &arguments->sorthdrnum, false);
  gt_option_exclude(sorthdrnum_option, invert_option);
  gt_option_exclude(sorthdrnum_option, sort_option);
  gt_option_exclude(sorthdrnum_option, revsort_option);
  gt_option_exclude(sorthdrnum_option, sorthdr_option);
  gt_option_parser_add_option(op, sorthdrnum_option);

  /* -shuffle */
  shuffle_option = gt_option_new_bool("shuffle", "shuffle sequences "
                                      "pseudo-randomly",
                                      &arguments->shuffle, false);
  gt_option_exclude(shuffle_option, invert_option);
  gt_option_exclude(shuffle_option, sort_option);
  gt_option_exclude(shuffle_option, revsort_option);
  gt_option_exclude(shuffle_option, sorthdr_option);
  gt_option_exclude(shuffle_option, sorthdrnum_option);
  gt_option_parser_add_option(op, shuffle_option);

  /* -sortlength */
  sort_length_option = gt_option_new_bool("sortlength",
                                          "sort by decreasing length",
                                          &arguments->sort_length, false);
  gt_option_exclude(sort_length_option, invert_option);
  gt_option_exclude(sort_length_option, sort_option);
  gt_option_exclude(sort_length_option, revsort_option);
  gt_option_exclude(sort_length_option, sorthdr_option);
  gt_option_exclude(sort_length_option, sorthdrnum_option);
  gt_option_exclude(sort_length_option, shuffle_option);
  gt_option_parser_add_option(op, sort_length_option);
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
        arguments->sorthdr || arguments->sorthdrnum || arguments->shuffle ||
        arguments->sort_length))
  {
    had_err = 1;
    gt_error_set(err, "order option needed: -invert|-sort|-revsort|"
                      "-sorthdr|-sorthdrnum|shuffle|sort_length");
  }
  return had_err;
}

static void gt_seqorder_sort(GtSuffixsortspace *suffixsortspace,
                             const GtEncseq *encseq)
{
  Sfxstrategy sfxstrategy;

  defaultsfxstrategy(&sfxstrategy, false);
  gt_suffixsortspace_init_seqstartpos(suffixsortspace,encseq);
  gt_sortallsuffixesfromstart(suffixsortspace,
                              gt_encseq_num_of_sequences(encseq),
                              encseq, GT_READMODE_FORWARD, NULL, 0,
                              &sfxstrategy, NULL, NULL, NULL);
}

static void gt_seqorder_get_shuffled_seqnums(GtUword nofseqs,
                                             GtUword *seqnums)
{
  GtUword i, j;

  gt_assert(seqnums != NULL);
  seqnums[0] = 0;
  for (i = 1UL; i < nofseqs; i++)
  {
    j = gt_rand_max(i);
    seqnums[i] = seqnums[j];
    seqnums[j] = i;
  }
}

static int seqorder_str_compare_lex(const void *v1, const void *v2, void *data)
{
  GtUword n1 = *(const GtUword*) v1,
          n2 = *(const GtUword*) v2,
          desclen1, desclen2;
  const char *desc1, *desc2;
  int rval;

  desc1 = gt_encseq_description((GtEncseq*) data, &desclen1, n1);
  desc2 = gt_encseq_description((GtEncseq*) data, &desclen2, n2);
  rval = strncmp(desc1, desc2, GT_MIN(desclen1, desclen2) * sizeof *desc1);
  if (rval == 0)
  {
    if (desclen1 > desclen2)
    {
      return 1;
    }
    if (desclen1 < desclen2)
    {
      return -1;
    }
    return 0;
  }
  return rval;
}

static int seqorder_str_compare_num(const void *v1, const void *v2, void *data)
{
  GtUword n1 = *(const GtUword*) v1,
          n2 = *(const GtUword*) v2,
          desclen1, desclen2,
          anum, bnum;
  int arval, brval, rval = 0;
  const char *desc1, *desc2;
  char buf[BUFSIZ];
  desc1 = gt_encseq_description((GtEncseq*) data, &desclen1, n1);
  desc2 = gt_encseq_description((GtEncseq*) data, &desclen2, n2);
  (void) strncpy(buf, desc1, GT_MIN(BUFSIZ, desclen1) * sizeof (char));
  buf[desclen1] = '\0';
  arval = gt_parse_uword(&anum, buf);
  (void) strncpy(buf, desc2, GT_MIN(BUFSIZ, desclen2) * sizeof (char));
  buf[desclen2] = '\0';
  brval = gt_parse_uword(&bnum, buf);
  if (arval == 0 && brval == 0)
    rval = anum-bnum;
  else if (arval == 0)
    return -1;
  else if (brval == 0)
    return 1;
  else
    rval = 0;
  return rval;
}

static int seqorder_length_compare(const void *v1, const void *v2, void *data)
{
  GtUword n1 = *(const GtUword*) v1,
          n2 = *(const GtUword*) v2;
  const GtUword *lengthtab = (const GtUword *) data;

  if (lengthtab[n1] < lengthtab[n2])
  {
    return 1;
  }
  if (lengthtab[n1] > lengthtab[n2])
  {
    return -1;
  }
  return 0;
}

static void gt_seqorder_output(GtUword seqnum,const GtEncseq *encseq)
{
  GtEncseqReader *esr;
  GtUword startpos, len, desclen = 0;
  const char *desc = NULL;
  GtUword i;

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

static void gt_seqorder_handle_cases(const GtEncseq *encseq,
                                     const GtSeqorderArguments *arguments)
{
  GtUword i, nofseqs = gt_encseq_num_of_sequences(encseq);
  GtSuffixsortspace *suffixsortspace = NULL;

  if (arguments->invert)
  {
    for (i = nofseqs; i > 0; i--)
    {
      gt_seqorder_output(i - 1, encseq);
    }
    return;
  }
  if (arguments->shuffle || arguments->sorthdr || arguments->sorthdrnum ||
      arguments->sort_length)
  {
    const GtUword numofsequences = gt_encseq_num_of_sequences(encseq);
    GtUword *seqnums = gt_malloc(sizeof (GtUword) * nofseqs);
    if (arguments->shuffle)
    {
      gt_seqorder_get_shuffled_seqnums(nofseqs, seqnums);
    } else
    {
      GtCompareWithData seqordercmpfunc;
      void *data;

      for (i = 0UL; i < numofsequences; i++)
      {
        seqnums[i] = i;
      }
      if (arguments->sort_length)
      {
        seqordercmpfunc = seqorder_length_compare;
        data = (void *) gt_all_sequence_lengths_get(encseq);
      } else
      {
        data = (void *) encseq;
        if (arguments->sorthdr)
        {
          seqordercmpfunc = seqorder_str_compare_lex;
        } else
        {
          if (arguments->sorthdrnum)
          {
            seqordercmpfunc = seqorder_str_compare_num;
          } else
          {
            gt_assert(false);
            seqordercmpfunc = NULL;
          }
        }
      }
      if (data != NULL)
      {
        gt_assert(seqordercmpfunc != NULL);
        (void) gt_qsort_r(seqnums, numofsequences,sizeof *seqnums,
                          data,seqordercmpfunc);
      }
      if (arguments->sort_length && data != NULL)
      {
        gt_free(data);
      }
    }
    for (i = 0; i < nofseqs; i++)
    {
      gt_seqorder_output(seqnums[i], encseq);
    }
    gt_free(seqnums);
    return;
  }
  gt_assert(arguments->sort || arguments->revsort);
  suffixsortspace
    = gt_suffixsortspace_new(nofseqs,
                             /* Use iterator over sequence separators:
                                saves a lot of binary searches */
                             gt_encseq_seqstartpos(encseq, nofseqs-1),
                             false,NULL);
  gt_seqorder_sort(suffixsortspace, encseq);
  if (arguments->sort)
  {
    for (i = 0; i < nofseqs; i++)
    {
      GtUword pos = gt_suffixsortspace_getdirect(suffixsortspace,i);
      gt_seqorder_output(gt_encseq_seqnum(encseq,pos),encseq);
    }
  } else
  {
    for (i = nofseqs; i > 0; i--)
    {
      GtUword pos = gt_suffixsortspace_getdirect(suffixsortspace,i - 1);
      gt_seqorder_output(gt_encseq_seqnum(encseq,pos),encseq);
    }
  }
  gt_suffixsortspace_delete(suffixsortspace, false);
}

static int gt_seqorder_runner(GT_UNUSED int argc,
                              const char **argv,
                              int parsed_args,
                              void *tool_arguments,
                              GtError *err)
{
  GtSeqorderArguments *arguments = tool_arguments;
  int had_err = 0;
  GtEncseq *encseq;
  GtEncseqLoader *loader;

  gt_error_check(err);
  gt_assert(arguments != NULL);

  /* load encseq */
  loader = gt_encseq_loader_new();
  encseq = gt_encseq_loader_load(loader, argv[parsed_args], err);
  gt_encseq_loader_delete(loader);
  if (encseq == NULL)
  {
    had_err = -1;
  }
  if (had_err == 0 && !gt_encseq_has_description_support(encseq))
  {
    gt_warning("%s has no description support", argv[parsed_args]);
  }
  if (!had_err)
  {
    gt_seqorder_handle_cases(encseq,arguments);
  }
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
