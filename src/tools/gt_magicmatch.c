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
#include "core/bioseq.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/unused_api.h"
#include "extended/gtdatahelp.h"
#include "tools/gt_magicmatch.h"

typedef struct {
  GtStrArray *seqfiles;
  bool translate;
} MagicMatchArguments;

static void* gt_magicmatch_arguments_new(void)
{
  MagicMatchArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->seqfiles = gt_str_array_new();
  return arguments;
}

static void gt_magicmatch_arguments_delete(void *tool_arguments)
{
  MagicMatchArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_array_delete(arguments->seqfiles);
  gt_free(arguments);
}

static GtOptionParser* gt_magicmatch_option_parser_new(void *tool_arguments)
{
  MagicMatchArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *o;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] -f sequence_file [...] -t",
                         "Compute MD5 fingerprints for each sequence given in "
                         "sequence_file(s).");

  /* -f */
  o = gt_option_new_filename_array("f", "fasta file names (at least one file "
                                   "is required)", arguments->seqfiles);
  gt_option_is_mandatory(o);
  gt_option_parser_add_option(op, o);

  /* -t */
  o = gt_option_new_bool("t", "translate the sequences of the files",
                         &arguments->translate, false);
  gt_option_is_mandatory(o);
  gt_option_parser_add_option(op, o);

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  gt_option_parser_set_min_max_args(op, 0, 0);

  return op;
}

static void translate_sequence_file(GtBioseq *bs)
{
  unsigned long i;
  gt_assert(bs);
  for (i = 0; i < gt_bioseq_number_of_sequences(bs); i++) {
    printf("%s\t%s\n", gt_bioseq_get_md5_fingerprint(bs, i),
                       gt_bioseq_get_description(bs, i));
  }
}

static int gt_magicmatch_runner(GT_UNUSED int argc, GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args, void *tool_arguments,
                                GtError *err)
{
  MagicMatchArguments *arguments = tool_arguments;
  GtBioseq *bioseq;
  unsigned long i;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (arguments->translate) {
    for (i = 0; !had_err && i < gt_str_array_size(arguments->seqfiles); i++) {
      if (!(bioseq = gt_bioseq_new(gt_str_array_get(arguments->seqfiles, i),
                                   err))) {
        had_err = -1;
      }
      if (!had_err)
        translate_sequence_file(bioseq);
      gt_bioseq_delete(bioseq);
    }
  }

  return had_err;
}

GtTool* gt_magicmatch(void)
{
  return gt_tool_new(gt_magicmatch_arguments_new,
                  gt_magicmatch_arguments_delete,
                  gt_magicmatch_option_parser_new,
                  NULL,
                  gt_magicmatch_runner);
}
