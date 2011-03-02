/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/parseutils_api.h"
#include "core/readmode_api.h"
#include "core/unused_api.h"
#include "core/undef.h"
#include "tools/gt_bitextract.h"

typedef struct {
  bool mirror,
       stoppos,
       specialranges;
  GtStr *readmode;
  unsigned long bitpos,
                iterations;
} GtBitextractArguments;

static void* gt_bitextract_arguments_new(void)
{
  GtBitextractArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->readmode = gt_str_new();
  return arguments;
}

static void gt_bitextract_arguments_delete(void *tool_arguments)
{
  GtBitextractArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->readmode);
  gt_free(arguments);
}

static GtOptionParser* gt_bitextract_option_parser_new(void *tool_arguments)
{
  GtBitextractArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] [indexname]",
                            "Extracts internal data from encoded sequences.");

  option = gt_option_new_bool("mirrored", "mirror sequence",
                               &arguments->mirror, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("stoppos", "output stop positions",
                               &arguments->stoppos, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("specialranges", "output special ranges",
                               &arguments->specialranges, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong("bitpos", "extract and display "
                                         "two bit encoding for position",
                               &arguments->bitpos, GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong("iterations", "number of test iterations",
                               &arguments->iterations, 10);
  gt_option_parser_add_option(op, option);

  gt_encseq_options_add_readmode_option(op, arguments->readmode);

  gt_option_parser_set_min_args(op, 1U);

  return op;
}

static int gt_bitextract_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GT_UNUSED GtBitextractArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
}

static int gt_bitextract_runner(GT_UNUSED int argc, const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments, GT_UNUSED GtError *err)
{
  GtBitextractArguments *arguments = tool_arguments;
  GtEncseqLoader *el;
  GtEncseq *encseq;
  int had_err = 0;
  bool fwd, it1, GT_UNUSED it2;
  char buffer[BUFSIZ];
  GtEndofTwobitencoding etbe;
  GtEncseqReader *esr;
  GtSpecialrangeiterator *sri;
  GtRange srng;
  GtReadmode rm;

  gt_error_check(err);
  gt_assert(arguments);

  el = gt_encseq_loader_new();
  encseq = gt_encseq_loader_load(el, argv[parsed_args], err);
  if (!encseq)
    had_err = -1;

  if (!had_err && arguments->mirror) {
    gt_encseq_mirror(encseq);
  }

  if (!had_err) {
    rm = gt_readmode_parse(gt_str_get(arguments->readmode), NULL);
    fwd = GT_ISDIRREVERSE(rm) ? false : true;
  }

  if (!had_err && arguments->bitpos != GT_UNDEF_ULONG) {
    if (arguments->bitpos >= gt_encseq_total_length(encseq)) {
      gt_error_set(err, "position %lu exceeds encoded sequence length of %lu",
                   arguments->bitpos, gt_encseq_total_length(encseq));
      had_err = -1;
    }

    if (!had_err) {
      unsigned long ret;
      esr = gt_encseq_create_reader_with_readmode(encseq, rm,
                                                  arguments->bitpos);
      ret = gt_encseq_extract2bitencwithtwobitencodingstoppos(&etbe, esr,
                                                        encseq,
                                                        rm, arguments->bitpos);
      gt_bitsequence_tostring(buffer, etbe.tbe);
      printf("Twobitencoding   %s\n"
             "unitsnotspecial  %u\n"
             "position         %lu\n"
             "returnvalue      %lu\n",
             buffer,
             etbe.unitsnotspecial,
             etbe.position,
             ret);
      gt_encseq_reader_delete(esr);
    }
  }

  if (!had_err && arguments->stoppos) {
    unsigned long i, stoppos;
    esr = gt_encseq_create_reader_with_readmode(encseq, rm, 0);
    /* check stoppos stuff */
    for (i=0;i<arguments->iterations;i++) {
      unsigned long pos;
      pos = gt_rand_max(gt_encseq_total_length(encseq) - 1);
      gt_encseq_reader_reinit_with_readmode(esr, encseq, rm, pos);
      stoppos = gt_getnexttwobitencodingstoppos(fwd, esr);
      printf("%lu: %lu\n", pos, stoppos);
    }
    gt_encseq_reader_delete(esr);
  }

  if (!had_err && arguments->specialranges) {
    /* check specialrangeiterator stuff */
    sri = gt_specialrangeiterator_new(encseq, fwd);
    while (true) {
      it1 = gt_specialrangeiterator_next(sri, &srng);
      if (it1)
        printf("%lu:%lu\n", srng.start, srng.end);
      else break;
    }
    gt_specialrangeiterator_delete(sri);
  }

  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  return had_err;
}

GtTool* gt_bitextract(void)
{
  return gt_tool_new(gt_bitextract_arguments_new,
                  gt_bitextract_arguments_delete,
                  gt_bitextract_option_parser_new,
                  gt_bitextract_arguments_check,
                  gt_bitextract_runner);
}
