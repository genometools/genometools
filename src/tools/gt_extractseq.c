/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/bioseq.h"
#include "libgtcore/fasta.h"
#include "libgtcore/grep.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/outputfile.h"
#include "tools/gt_extractseq.h"

typedef struct {
  Str *pattern;
  unsigned long width;
  OutputFileInfo *ofi;
  GenFile *outfp;
} ExtractSeqArguments;

static void* gt_extractseq_arguments_new(void)
{
  ExtractSeqArguments *arguments = ma_calloc(1, sizeof *arguments);
  arguments->pattern = str_new();
  arguments->ofi = outputfileinfo_new();
  return arguments;
}

static void gt_extractseq_arguments_delete(void *tool_arguments)
{
  ExtractSeqArguments *arguments = tool_arguments;
  if (!tool_arguments) return;
  genfile_close(arguments->outfp);
  outputfileinfo_delete(arguments->ofi);
  str_delete(arguments->pattern);
  ma_free(arguments);
}

static OptionParser* gt_extractseq_option_parser_new(void *tool_arguments)
{
  ExtractSeqArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *option;
  assert(arguments);

  /* init */
  op = option_parser_new("[option ...] [sequence_file ...]",
                         "Extract sequences from given sequence file(s).");

  /* -match */
  option = option_new_string("match", "extract all sequences whose description "
                             "matches the given pattern.\nThe given pattern "
                             "must be a valid extended regular expression.",
                             arguments->pattern, NULL);
  option_is_mandatory(option);
  option_parser_add_option(op, option);

  /* -width */
  option = option_new_ulong("width", "set output width for showing of "
                            "sequences (0 disables formatting)",
                            &arguments->width, 0);
  option_parser_add_option(op, option);

  /* output file options */
  outputfile_register_options(op, &arguments->outfp, arguments->ofi);

  return op;
}

static int extractseq(GenFile *outfp, Bioseq *bs, const char *pattern,
                      unsigned long width, Error *err)
{
  const char *desc;
  unsigned long i;
  bool match;
  int had_err = 0;

  error_check(err);
  assert(bs && pattern);

  for (i = 0; !had_err && i < bioseq_number_of_sequences(bs); i++) {
    desc = bioseq_get_description(bs, i);
    assert(desc);
    had_err = grep(&match, pattern, desc, err);
    if (!had_err && match) {
      fasta_show_entry_generic(desc, bioseq_get_sequence(bs, i),
                               bioseq_get_sequence_length(bs, i), width, outfp);
    }
  }

  return had_err;
}

static int gt_extractseq_runner(int argc, const char **argv,
                                void *tool_arguments, Error *err)
{
  ExtractSeqArguments *arguments = tool_arguments;
  Bioseq *bs;
  int arg = 0, had_err = 0;

  error_check(err);
  assert(arguments);

  if (argc == 0) { /* no file given, use stdin */
    if (!(bs = bioseq_new("-", err)))
      had_err = -1;
    if (!had_err) {
      had_err = extractseq(arguments->outfp, bs, str_get(arguments->pattern),
                           arguments->width, err);
    }
    bioseq_delete(bs);
  }

  /* process all files */
  while (!had_err && arg < argc) {
    if (!(bs = bioseq_new(argv[arg], err)))
      had_err = -1;
    if (!had_err) {
      had_err = extractseq(arguments->outfp, bs, str_get(arguments->pattern),
                           arguments->width, err);
    }
    bioseq_delete(bs);
    arg++;
  }

  return had_err;
}

Tool* gt_extractseq(void)
{
  return tool_new(gt_extractseq_arguments_new,
                  gt_extractseq_arguments_delete,
                  gt_extractseq_option_parser_new,
                  NULL,
                  gt_extractseq_runner);
}
