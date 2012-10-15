/*
  Copyright (c) 2010      Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2010-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2012 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/ma.h"
#include "core/chardef.h"
#include "core/encseq_api.h"
#include "core/encseq_options.h"
#include "core/fasta_separator.h"
#include "core/log_api.h"
#include "core/readmode.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "tools/gt_encseq_decode.h"

typedef struct {
  bool singlechars;
  GtStr *mode,
        *sepchar;
  GtRange rng,
          seqrng;
  GtEncseqOptions *eopts;
  GtReadmode rm;
  GtStr *dir;
  unsigned long seq;
} GtEncseqDecodeArguments;

static void* gt_encseq_decode_arguments_new(void)
{
  GtEncseqDecodeArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->mode = gt_str_new();
  arguments->sepchar = gt_str_new();
  arguments->dir = gt_str_new();
  arguments->seqrng.start = arguments->seqrng.end = GT_UNDEF_ULONG;
  arguments->rng.start = arguments->rng.end = GT_UNDEF_ULONG;
  return arguments;
}

static void gt_encseq_decode_arguments_delete(void *tool_arguments)
{
  GtEncseqDecodeArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->mode);
  gt_str_delete(arguments->sepchar);
  gt_str_delete(arguments->dir);
  gt_encseq_options_delete(arguments->eopts);
  gt_free(arguments);
}

static GtOptionParser* gt_encseq_decode_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option,
           *optionsep,
           *optionseq,
           *optionseqrange,
           *optionmode;
  GtEncseqDecodeArguments *arguments =
                                      (GtEncseqDecodeArguments*) tool_arguments;
  static const char *modes[] = {"fasta", "concat", NULL};

  /* init */
  op = gt_option_parser_new("(sequence_file|indexname)",
                            "Decode/extract encoded sequences.");

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

  /* -seq */
  optionseq = gt_option_new_ulong("seq",
                                  "extract sequence identified by its number",
                                  &arguments->seq,
                                  GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, optionseq);

  /* -seqrange */
  optionseqrange = gt_option_new_range("seqrange",
                                       "extract multiple consecutive sequences",
                                       &arguments->seqrng,
                                       NULL);
  gt_option_parser_add_option(op, optionseqrange);
  gt_option_exclude(optionseqrange, optionseq);
  gt_option_exclude(optionseq, optionseqrange);

  /* -output */
  optionmode = gt_option_new_choice("output",
                                    "specify output format "
                                    "(choose from fasta|concat)",
                                    arguments->mode,
                                    modes[0],
                                    modes);
  gt_option_parser_add_option(op, optionmode);

  /* -range */
  option = gt_option_new_range("range",
                               "concatenated range to extract "
                               "(implies '-output concat')",
                               &arguments->rng,
                               NULL);
  gt_option_parser_add_option(op, option);
  gt_option_imply(option, optionmode);

  /* -sepchar */
  optionsep = gt_option_new_string("sepchar",
                                   "specify character to print as SEPARATOR",
                                   arguments->sepchar, "|");
  gt_option_parser_add_option(op, optionsep);
  gt_option_imply(optionsep, optionmode);

  gt_option_parser_set_min_max_args(op, 1, 1);

  return op;
}

int gt_encseq_decode_arguments_check(GT_UNUSED int rest_argc,
                                     void *tool_arguments,
                                     GtError* err)
{
  GtEncseqDecodeArguments *args = (GtEncseqDecodeArguments*) tool_arguments;
  int had_err = 0;
  int rval;

  if (gt_str_length(args->dir) > 0) {
    rval = gt_readmode_parse(gt_str_get(args->dir), err);
    if (rval < 0)
      had_err = -1;
    else
      args->rm = (GtReadmode) rval;
  }
  if (!had_err && args->seqrng.start != GT_UNDEF_ULONG &&
        args->seqrng.end != GT_UNDEF_ULONG && args->seq != GT_UNDEF_ULONG) {
    gt_error_set(err, "'-seq' and '-seqrange' cannot be used together");
    had_err = -1;
  }
  if (!had_err && (args->seqrng.start != GT_UNDEF_ULONG ||
        args->seqrng.end != GT_UNDEF_ULONG || args->seq != GT_UNDEF_ULONG)
        && strcmp(gt_str_get(args->mode), "fasta") != 0) {
    gt_error_set(err, "'-seq' and '-seqrange' can only be used with the "
                      "'-output fasta' option");
    had_err = -1;
  }
  if (!had_err && (args->rng.start != GT_UNDEF_ULONG ||
        args->rng.end != GT_UNDEF_ULONG)
        && strcmp(gt_str_get(args->mode), "concat") != 0) {
    gt_error_set(err, "'-range' can only be used with the "
                      "'-output concat' option");
    had_err = -1;
  }
  return had_err;
}

static int output_sequence(GtEncseq *encseq, GtEncseqDecodeArguments *args,
                           const char *filename, GtError *err)
{
  unsigned long i, j, sfrom, sto;
  int had_err = 0;
  bool has_desc;
  GtEncseqReader *esr;
  gt_assert(encseq);

  if (!(has_desc = gt_encseq_has_description_support(encseq)))
    gt_warning("Missing description support for file %s", filename);

  if (strcmp(gt_str_get(args->mode), "fasta") == 0) {
    /* specify a single sequence to extract */
    if (args->seq != GT_UNDEF_ULONG) {
      if (args->seq >= gt_encseq_num_of_sequences(encseq)) {
        gt_error_set(err, "requested sequence %lu exceeds number of sequences "
                          "(%lu)", args->seq,
                          gt_encseq_num_of_sequences(encseq));
        return -1;
      }
      sfrom = args->seq;
      sto = args->seq + 1;
    } else if (args->seqrng.start != GT_UNDEF_ULONG
                 && args->seqrng.end != GT_UNDEF_ULONG) {
      /* specify a sequence range to extract */
      if (args->seqrng.start >= gt_encseq_num_of_sequences(encseq)
            || args->seqrng.end >= gt_encseq_num_of_sequences(encseq)) {
        gt_error_set(err, "range %lu-%lu includes a sequence number exceeding "
                          "the total number of sequences (%lu)",
                          args->seqrng.start,
                          args->seqrng.end,
                          gt_encseq_num_of_sequences(encseq));
        return -1;
      }
      sfrom = args->seqrng.start;
      sto = args->seqrng.end + 1;
    } else {
      /* extract all sequences */
      sfrom = 0;
      sto = gt_encseq_num_of_sequences(encseq);
    }
    for (i = sfrom; i < sto; i++) {
      unsigned long desclen, startpos, len;
      char buf[BUFSIZ];
      const char *desc = NULL;
      /* XXX: maybe make this distinction in the functions via readmode? */
      if (!GT_ISDIRREVERSE(args->rm)) {
        startpos = gt_encseq_seqstartpos(encseq, i);
        len = gt_encseq_seqlength(encseq, i);
        if (has_desc) {
          desc = gt_encseq_description(encseq, &desclen, i);
        } else {
          (void) snprintf(buf, BUFSIZ, "sequence %lu", i);
          desclen = strlen(buf);
          desc = buf;
        }
      } else {
        startpos = gt_encseq_seqstartpos(encseq, i);
        len = gt_encseq_seqlength(encseq,
                                  gt_encseq_num_of_sequences(encseq)-1-i);
        startpos = gt_encseq_total_length(encseq)
                     - (gt_encseq_seqstartpos(encseq,
                                              gt_encseq_num_of_sequences(
                                                encseq)-1-i) + len);
        if (has_desc) {
          desc = gt_encseq_description(encseq,
                                       &desclen,
                                       gt_encseq_num_of_sequences(encseq)-1-i);
        } else {
          (void) snprintf(buf, BUFSIZ, "sequence %lu", i);
          desclen = strlen(buf);
          desc = buf;
        }
      }
      gt_assert(desc);
      /* output description */
      gt_xfputc(GT_FASTA_SEPARATOR, stdout);
      gt_xfwrite(desc, 1, desclen, stdout);
      gt_xfputc('\n', stdout);
      /* XXX: make this more efficient by writing in a buffer first and then
         showing the result */
      if (args->singlechars) {
        for (j = 0; j < len; j++) {
           gt_xfputc(gt_encseq_get_decoded_char(encseq,
                                                startpos + j,
                                                args->rm),
                     stdout);
        }
      } else {
        esr = gt_encseq_create_reader_with_readmode(encseq, args->rm, startpos);
        for (j = 0; j < len; j++) {
           gt_xfputc(gt_encseq_reader_next_decoded_char(esr), stdout);
        }
        gt_encseq_reader_delete(esr);
      }
      gt_xfputc('\n', stdout);
    }
  }

  if (strcmp(gt_str_get(args->mode), "concat") == 0) {
    unsigned long from = 0,
                  to = gt_encseq_total_length(encseq) - 1;
    if (args->rng.start != GT_UNDEF_ULONG && args->rng.end != GT_UNDEF_ULONG) {
      if (args->rng.end > to) {
        had_err = -1;
        gt_error_set(err, "end of range (%lu) exceeds encoded sequence length "
                          "(%lu)", args->rng.end, to);
      }
      if (!had_err) {
        from = args->rng.start;
        to = args->rng.end;
      }
    }
    if (!had_err) {
      if (args->singlechars) {
        for (j = from; j <= to; j++) {
          char cc = gt_encseq_get_decoded_char(encseq, j, args->rm);
          if (cc == (char) SEPARATOR)
            cc = gt_str_get(args->sepchar)[0];
          gt_xfputc(cc, stdout);
        }
      } else {
        esr = gt_encseq_create_reader_with_readmode(encseq, args->rm, from);
        if (esr) {
          for (j = from; j <= to; j++) {
            char cc = gt_encseq_reader_next_decoded_char(esr);
            if (cc == (char) SEPARATOR)
              cc = gt_str_get(args->sepchar)[0];
            gt_xfputc(cc, stdout);
          }
          gt_encseq_reader_delete(esr);
        }
      }
      gt_xfputc('\n', stdout);
    }
  }
  return had_err;
}

static int decode_sequence_file(const char *seqfile,
                                GtEncseqDecodeArguments *args,
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
  if (!had_err)
    had_err = output_sequence(encseq, args, seqfile, err);
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(encseq_loader);
  return had_err;
}

static int gt_encseq_decode_runner(GT_UNUSED int argc, const char **argv,
                               int parsed_args, void *tool_arguments,
                               GtError *err)
{
  gt_error_check(err);
  return decode_sequence_file(argv[parsed_args], tool_arguments, err);
}

GtTool* gt_encseq_decode(void)
{
  return gt_tool_new(gt_encseq_decode_arguments_new,
                     gt_encseq_decode_arguments_delete,
                     gt_encseq_decode_option_parser_new,
                     gt_encseq_decode_arguments_check,
                     gt_encseq_decode_runner);
}
