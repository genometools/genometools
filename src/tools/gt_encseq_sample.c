/*
  Copyright (c) 2014 Joerg Winkler <j.winkler@posteo.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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
#include "tools/gt_encseq_sample.h"

#include "core/chardef.h"
#include "core/encseq_api.h"
#include "core/encseq_options.h"
#include "core/fasta_separator.h"
#include "core/log_api.h"
#include "core/readmode.h"
#include "core/undef_api.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "core/minmax.h"
#include "core/mathsupport.h"
#include "core/bittab_api.h"

#include <string.h>
#include <stdlib.h>

typedef struct {
  bool singlechars;
  GtStr *mode,
        *sepchar;
  GtRange seqrng;
  GtEncseqOptions *eopts;
  GtReadmode rm;
  GtStr *dir;
  GtUword len;
} GtEncseqSampleArguments;

typedef struct {
  GtUword offset;
  GtBittab *sample;
} GtEncseqSampleOutputInfo;

static void* gt_encseq_sample_arguments_new(void)
{
  GtEncseqSampleArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->mode = gt_str_new();
  arguments->sepchar = gt_str_new();
  arguments->dir = gt_str_new();
  arguments->seqrng.start = arguments->seqrng.end = GT_UNDEF_UWORD;
  return arguments;
}

static void gt_encseq_sample_arguments_delete(void *tool_arguments)
{
  GtEncseqSampleArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->mode);
    gt_str_delete(arguments->sepchar);
    gt_str_delete(arguments->dir);
    gt_encseq_options_delete(arguments->eopts);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_encseq_sample_option_parser_new(void *tool_arguments)
{
  GtEncseqSampleArguments *arguments = (GtEncseqSampleArguments*)tool_arguments;
  GtOptionParser *op;
  GtOption *option,
           *optionsep,
           *optionlen,
           *optionseqrange,
           *optionmode;
  gt_assert(arguments);
  static const char *modes[] = {"fasta", "concat", NULL};

  /* init */
  op = gt_option_parser_new("(sequence_file|indexname)",
                            "Decode/extract encoded sequences by random "
                            "choice.");

  /* encseq options */
  arguments->eopts = gt_encseq_options_register_loading(op, NULL);
  gt_encseq_options_add_readmode_option(op, arguments->dir);

  /* -singlechars */
  option = gt_option_new_bool("singlechars",
                              "do not use a GtEncseqReader but access each "
                              "sequence character separately",
                              &arguments->singlechars,
                              false);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* -length */
  optionlen = gt_option_new_uword("length",
                                  "minimum length to be extracted",
                                  &arguments->len,
                                  GT_UNDEF_UWORD);
  gt_option_parser_add_option(op, optionlen);

  /* -seqrange */
  optionseqrange = gt_option_new_range("seqrange",
                                       "extract multiple consecutive sequences",
                                       &arguments->seqrng,
                                       NULL);
  gt_option_parser_add_option(op, optionseqrange);

  /* -output */
  optionmode = gt_option_new_choice("output",
                                    "specify output format "
                                    "(choose from fasta|concat)",
                                    arguments->mode,
                                    modes[0],
                                    modes);
  gt_option_parser_add_option(op, optionmode);

  /* -sepchar */
  optionsep = gt_option_new_string("sepchar",
                                   "specify character to print as SEPARATOR",
                                   arguments->sepchar, "|");
  gt_option_parser_add_option(op, optionsep);
  gt_option_imply(optionsep, optionmode);

  gt_option_parser_set_min_max_args(op, 1, 1);

  return op;
}

static int gt_encseq_sample_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtEncseqSampleArguments *args = (GtEncseqSampleArguments *) tool_arguments;
  int rval;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(args);

  if (gt_str_length(args->dir) > 0) {
    rval = gt_readmode_parse(gt_str_get(args->dir), err);
    if (rval < 0)
      had_err = -1;
    else
      args->rm = (GtReadmode) rval;
  }
  if (args->len == GT_UNDEF_UWORD) {
    gt_error_set(err, "length must be specified");
    return -1;
  }
  return had_err;
}

static int gt_encseq_sample_output(GtEncseq *encseq,
                                   GtEncseqSampleArguments *args,
                                   GtEncseqSampleOutputInfo *output)
{
  GtUword i, j, stop;
  GtEncseqReader *esr;
  const bool is_concat = (strcmp(gt_str_get(args->mode), "concat") == 0);
  const bool is_reverse = GT_ISDIRREVERSE(args->rm);
  int had_err = 0;

  /* extract selected sequences */
  i = gt_bittab_get_first_bitnum(output->sample) + output->offset;
  stop = gt_bittab_get_last_bitnum(output->sample) + output->offset;
  while (i < stop) {
    GtUword startpos, len;

    if (is_reverse) {
      len = gt_encseq_seqlength(encseq,
                                gt_encseq_num_of_sequences(encseq) - 1 - i);
      startpos = gt_encseq_total_length(encseq) - (gt_encseq_seqstartpos(encseq,
                 gt_encseq_num_of_sequences(encseq) - 1 - i) + len);
    } else {
      startpos = gt_encseq_seqstartpos(encseq, i);
      len = gt_encseq_seqlength(encseq, i);
    }
    /* prepare description */
    if (!is_concat) {
      GtUword desclen;
      const char *desc = NULL;
      if (gt_encseq_has_description_support(encseq)) {
        if (is_reverse) {
          desc = gt_encseq_description(encseq, &desclen,
                                       gt_encseq_num_of_sequences(encseq)-1-i);
        } else {
          desc = gt_encseq_description(encseq, &desclen, i);
        }
      } else {
        char buf[BUFSIZ];
        (void) snprintf(buf, BUFSIZ, "sequence "GT_WU"", i);
        desclen = strlen(buf);
        desc = buf;
      }
      gt_assert(desc);
      /* output description */
      gt_xfputc(GT_FASTA_SEPARATOR, stdout);
      gt_xfwrite(desc, 1, desclen, stdout);
      gt_xfputc('\n', stdout);
    }

    if (args->singlechars) {
      for (j = 0; j < len; j++) {
        gt_xfputc(gt_encseq_get_decoded_char(encseq, startpos + j, args->rm),
                  stdout);
      }
    } else {
      esr = gt_encseq_create_reader_with_readmode(encseq, args->rm, startpos);
      for (j = 0; j < len; j++) {
        gt_xfputc(gt_encseq_reader_next_decoded_char(esr), stdout);
      }
      gt_encseq_reader_delete(esr);
    }
    i = gt_bittab_get_next_bitnum(output->sample, i - output->offset)
        + output->offset;

    /* insert separator between sequences */
    if (is_concat && i < stop) {
      gt_xfputc(gt_str_get(args->sepchar)[0], stdout);
    } else {
      gt_xfputc('\n', stdout);
    }
  }
  gt_bittab_delete(output->sample);
  return had_err;
}

static int gt_encseq_sample_choose_sequences(GtEncseq *encseq,
                                             GtEncseqSampleArguments *args,
                                             GtError *err,
                                             GtEncseqSampleOutputInfo *output)
{
  GtUword i, sfrom, sto;
  GtUword num_sequences, count;
  GtUword sequence_length;
  GtUword total_num_seq;
  int had_err = 0;

  gt_assert(encseq);
  sequence_length = gt_encseq_min_seq_length(encseq);
  total_num_seq = gt_encseq_num_of_sequences(encseq);

  if (sequence_length != gt_encseq_max_seq_length(encseq)) {
    gt_error_set(err, "sequences do not have the same length");
    return -1;
  }
  if (args->seqrng.start != GT_UNDEF_UWORD &&
      args->seqrng.end != GT_UNDEF_UWORD) {
    /* specify a sequence range to extract */
    if (args->seqrng.start > args->seqrng.end) {
      gt_error_set(err, "range start ("GT_WU") must not be higher than "
      "range end ("GT_WU")",
      args->seqrng.start, args->seqrng.end);
      return -1;
    }
    if (args->seqrng.end >= total_num_seq) {
      gt_error_set(err, "range "GT_WU"-"GT_WU" includes a sequence number "
                   "exceeding the total number of sequences ("GT_WU")",
                   args->seqrng.start, args->seqrng.end,
                   gt_encseq_num_of_sequences(encseq));
      return -1;
    }
    sfrom = args->seqrng.start;
    sto = args->seqrng.end;
    total_num_seq = 1+sto-sfrom;
  } else {
    /* extract all sequences */
    sfrom = 0;
    sto = total_num_seq-1;
  }
  output->offset = sfrom;
  if (args->len > total_num_seq * sequence_length) {
    gt_error_set(err, "requested length "GT_WU" exceeds length of sequences"
    " ("GT_WU")", args->len, total_num_seq*sequence_length);
    return -1;
  }

  /* fill bit vector randomly */
  output->sample = gt_bittab_new(total_num_seq);
  num_sequences = ceil(args->len / (double)sequence_length);
  if (total_num_seq != 1) {
    count = 0;
    i = gt_rand_max(total_num_seq-1);
    while (count < num_sequences) {
      /* probability = num_sequences/total_num_seq */
      if (gt_rand_max(total_num_seq-1) < num_sequences &&
        !gt_bittab_bit_is_set(output->sample, i)) {
        gt_bittab_set_bit(output->sample, i);
      count ++;
        }
        i = (i+1) % total_num_seq;
    }
  } else {
    gt_bittab_set_bit(output->sample, 0);
  }
  return had_err;
}

static GtEncseq *gt_encseq_sample_get_encseq(const char *seqfile,
                                             GtEncseqSampleArguments *args,
                                             GtError *err)
{
  GtEncseqLoader *encseq_loader;
  GtEncseq *encseq;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(seqfile);
  encseq_loader = gt_encseq_loader_new();
  if (!had_err && gt_encseq_options_lossless_value(args->eopts)) {
    gt_encseq_loader_require_lossless_support(encseq_loader);
  }
  if (!(encseq = gt_encseq_loader_load(encseq_loader, seqfile, err)))
    had_err = -1;
  if (!had_err && gt_encseq_options_mirrored_value(args->eopts)) {
    if (!gt_alphabet_is_dna(gt_encseq_alphabet(encseq))) {
      gt_error_set(err, "mirroring is only defined on DNA sequences");
      had_err = -1;
    }
    if (!had_err)
      had_err = gt_encseq_mirror(encseq, err);
  }
  gt_encseq_loader_delete(encseq_loader);
  if (!had_err) {
    if (!gt_encseq_has_description_support(encseq))
      gt_warning("Missing description support for file %s", seqfile);
    return encseq;
  } else {
    return NULL;
  }
}

static int gt_encseq_sample_runner(GT_UNUSED int argc, const char **argv,
                                   int parsed_args, void *tool_arguments,
                                   GT_UNUSED GtError *err)
{
  GtEncseq *encseq;
  GtEncseqSampleOutputInfo output;
  int had_err = 0;
  gt_error_check(err);
  encseq = gt_encseq_sample_get_encseq(argv[parsed_args], tool_arguments, err);
  if (encseq != NULL)
    had_err = gt_encseq_sample_choose_sequences(encseq, tool_arguments, err,
                                                &output);
  else
    had_err = -1;
  if (!had_err)
    had_err = gt_encseq_sample_output(encseq, tool_arguments, &output);
  gt_encseq_delete(encseq);
  return had_err;
}

GtTool* gt_encseq_sample(void)
{
  return gt_tool_new(gt_encseq_sample_arguments_new,
                     gt_encseq_sample_arguments_delete,
                     gt_encseq_sample_option_parser_new,
                     gt_encseq_sample_arguments_check,
                     gt_encseq_sample_runner);
}
