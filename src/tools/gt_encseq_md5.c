/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include "core/encseq.h"
#include "core/ma.h"
#include "core/md5_fingerprint_api.h"
#include "core/output_file_api.h"
#include "core/unused_api.h"
#include "tools/gt_encseq_md5.h"

typedef struct {
  GtOutputFileInfo *ofi;
  GtFile *outfp;
  bool fromindex;
} GtEncseqInfoArguments;

static void* gt_encseq_md5_arguments_new(void)
{
  GtEncseqInfoArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_encseq_md5_arguments_delete(void *tool_arguments)
{
  GtEncseqInfoArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_encseq_md5_option_parser_new(void *tool_arguments)
{
  GtEncseqInfoArguments *arguments = tool_arguments;
  GtOption *option;
  GtOptionParser *op;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] indexname",
                            "Display MD5 sums for an encoded sequence.");

  /* -fromindex */
  option = gt_option_new_bool("fromindex", "use MD5 table from .md5 file",
                              &arguments->fromindex, true);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_min_max_args(op, 1, 1);
  return op;
}

static int gt_encseq_md5_runner(GT_UNUSED int argc, const char **argv,
                           int parsed_args, void *tool_arguments,
                           GtError *err)
{
  GtEncseqInfoArguments *arguments = tool_arguments;
  int had_err = 0;
  GtEncseqLoader *encseq_loader;
  GtEncseq *encseq;
  gt_error_check(err);
  gt_assert(arguments);

  encseq_loader = gt_encseq_loader_new();
  if (!(encseq = gt_encseq_loader_load(encseq_loader,
                                       argv[parsed_args], err)))
    had_err = -1;

  if (!had_err) {
    unsigned long i;

    if (arguments->fromindex) {
      GtMD5Tab *tab;

      if (!gt_encseq_has_md5_support(encseq)) {
        gt_error_set(err, "encoded sequence given by indexname \"%s\" does not "
                          "have MD5 support", argv[parsed_args]);
        had_err = -1;
      }
      if (!had_err) {
        tab = gt_encseq_get_md5_tab(encseq, err);
        if (tab) {
          for (i = 0; i < gt_encseq_num_of_sequences(encseq); i++) {
            gt_file_xprintf(arguments->outfp,
                            "%lu: %s\n", i, gt_md5_tab_get(tab, i));
          }
          gt_md5_tab_delete(tab);
        } else had_err = -1;
      }
    } else {
      char *seq,
           *md5str;
      for (i = 0; i < gt_encseq_num_of_sequences(encseq); i++) {
        unsigned long len, start, end;
        len = gt_encseq_seqlength(encseq, i);
        start = gt_encseq_seqstartpos(encseq, i);
        end = len + start -1;
        seq = gt_malloc(len * sizeof (char));
        gt_encseq_extract_decoded(encseq, seq, start, end);
        md5str = gt_md5_fingerprint(seq, len);
        gt_file_xprintf(arguments->outfp, "%lu: %s\n", i, md5str);
        gt_free(seq);
        gt_free(md5str);
      }
    }
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(encseq_loader);
  return had_err;
}

GtTool* gt_encseq_md5(void)
{
  return gt_tool_new(gt_encseq_md5_arguments_new,
                  gt_encseq_md5_arguments_delete,
                  gt_encseq_md5_option_parser_new,
                  NULL,
                  gt_encseq_md5_runner);
}
