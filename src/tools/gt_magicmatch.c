/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/bioseq.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/unused.h"
#include "libgtext/gtdatahelp.h"
#include "tools/gt_magicmatch.h"

typedef struct {
  StrArray *seqfiles;
  bool translate;
} MagicMatchArguments;

static void* gt_magicmatch_arguments_new(void)
{
  MagicMatchArguments *arguments = ma_calloc(1, sizeof *arguments);
  arguments->seqfiles = strarray_new();
  return arguments;
}

static void gt_magicmatch_arguments_delete(void *tool_arguments)
{
  MagicMatchArguments *arguments = tool_arguments;
  if (!arguments) return;
  strarray_delete(arguments->seqfiles);
  ma_free(arguments);
}

static OptionParser* gt_magicmatch_option_parser_new(void *tool_arguments)
{
  MagicMatchArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *o;
  assert(arguments);

  /* init */
  op = option_parser_new("[option ...] -f sequence_file [...] -t",
                         "Compute MD5 fingerprints for each sequence given in "
                         "sequence_file(s).");

  /* -f */
  o = option_new_filenamearray("f", "fasta file names (at least one file is "
                               "required)", arguments->seqfiles);
  option_is_mandatory(o);
  option_parser_add_option(op, o);

  /* -t */
  o = option_new_bool("t", "translate the sequences of the files",
                      &arguments->translate, false);
  option_is_mandatory(o);
  option_parser_add_option(op, o);

  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  option_parser_set_min_max_args(op, 0, 0);

  return op;
}

static void translate_sequence_file(Bioseq *bs)
{
  unsigned long i;
  assert(bs);
  for (i = 0; i < bioseq_number_of_sequences(bs); i++) {
    printf("%s\t%s\n", bioseq_get_md5_fingerprint(bs, i),
                       bioseq_get_description(bs, i));
  }
}

static int gt_magicmatch_runner(UNUSED int argc, UNUSED const char **argv,
                                UNUSED int parsed_args, void *tool_arguments,
                                Error *err)
{
  MagicMatchArguments *arguments = tool_arguments;
  Bioseq *bioseq;
  unsigned long i;
  int had_err = 0;

  error_check(err);
  assert(arguments);

  if (arguments->translate) {
    for (i = 0; !had_err && i < strarray_size(arguments->seqfiles); i++) {
      if (!(bioseq = bioseq_new(strarray_get(arguments->seqfiles, i), err)))
        had_err = -1;
      if (!had_err)
        translate_sequence_file(bioseq);
      bioseq_delete(bioseq);
    }
  }

  return had_err;
}

Tool* gt_magicmatch(void)
{
  return tool_new(gt_magicmatch_arguments_new,
                  gt_magicmatch_arguments_delete,
                  gt_magicmatch_option_parser_new,
                  NULL,
                  gt_magicmatch_runner);
}
