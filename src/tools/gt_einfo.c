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
#include "core/outputfile.h"
#include "core/unused_api.h"
#include "tools/gt_einfo.h"

typedef struct {
  bool bool_option_einfo;
  GtStr  *str_option_einfo;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} GtEinfoArguments;

static void* gt_einfo_arguments_new(void)
{
  GtEinfoArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->str_option_einfo = gt_str_new();
  arguments->ofi = gt_outputfileinfo_new();
  return arguments;
}

static void gt_einfo_arguments_delete(void *tool_arguments)
{
  GtEinfoArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->str_option_einfo);
  gt_file_delete(arguments->outfp);
  gt_outputfileinfo_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_einfo_option_parser_new(void *tool_arguments)
{
  GtEinfoArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GT_UNUSED GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] indexname",
                            "Display meta-information about an "
                            "encoded sequence.");

  /* output file options */
  gt_outputfile_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_min_args(op, 1);
  return op;
}

static int gt_einfo_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtEinfoArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
}

static int gt_einfo_runner(GT_UNUSED int argc, const char **argv,
                           int parsed_args, void *tool_arguments,
                           GtError *err)
{
  GtEinfoArguments *arguments = tool_arguments;
  int had_err = 0;
  GtEncseqLoader *encseq_loader;
  GtEncseq *encseq;
  gt_error_check(err);
  gt_assert(arguments);

  encseq_loader = gt_encseq_loader_new();
  if (!(encseq = gt_encseq_loader_load(encseq_loader, argv[parsed_args], err)))
    had_err = -1;
  if (!had_err) {
    GtAlphabet *alpha;
    const GtUchar *chars;
    const GtStrArray *filenames;
    unsigned long i;

    gt_file_xprintf(arguments->outfp, "index name: ");
    gt_file_xprintf(arguments->outfp, "%s\n", argv[parsed_args]);

    gt_file_xprintf(arguments->outfp, "total length: ");
    gt_file_xprintf(arguments->outfp, "%lu\n", gt_encseq_total_length(encseq));

    gt_file_xprintf(arguments->outfp, "compressed size: ");
    gt_file_xprintf(arguments->outfp, "%lu bytes\n",
                                      gt_encseq_sizeofrep(encseq));

    gt_file_xprintf(arguments->outfp, "number of sequences: ");
    gt_file_xprintf(arguments->outfp, "%lu\n",
                                      gt_encseq_num_of_sequences(encseq));

    gt_file_xprintf(arguments->outfp, "number of files: ");
    gt_file_xprintf(arguments->outfp, "%lu\n", gt_encseq_num_of_files(encseq));

    filenames = gt_encseq_filenames(encseq);
    gt_file_xprintf(arguments->outfp, "original filenames:\n");
    for (i = 0; i < gt_str_array_size(filenames); i++) {
      gt_file_xprintf(arguments->outfp, "\t%s (%lu characters)\n",
                                        gt_str_array_get(filenames, i),
                                        (unsigned long)
                                     gt_encseq_effective_filelength(encseq, i));
    }

    alpha = gt_encseq_alphabet(encseq);
    chars = gt_alphabet_characters(alpha);
    gt_file_xprintf(arguments->outfp, "alphabet size: ");
    gt_file_xprintf(arguments->outfp, "%u\n", gt_alphabet_num_of_chars(alpha));
    gt_file_xprintf(arguments->outfp, "alphabet characters: ");
    gt_file_xprintf(arguments->outfp, "%.*s\n", gt_alphabet_num_of_chars(alpha),
                                      (char*) chars);

    gt_file_xprintf(arguments->outfp, "character distribution:\n");
    for (i = 0; i < gt_alphabet_num_of_chars(alpha); i++) {
      unsigned long cc;
      cc = gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, chars[i]));
      gt_file_xprintf(arguments->outfp, "\t%c: %lu (%.2f%%)\n",
                                        (char) chars[i],
                                        cc,
                             (cc /(double) gt_encseq_total_length(encseq))*100);
    }

    gt_file_xprintf(arguments->outfp, "number of wildcards: ");
    gt_file_xprintf(arguments->outfp, "%lu (%lu range(s))\n",
                                      gt_encseq_wildcards(encseq),
                                      gt_encseq_realwildcardranges(encseq));

        gt_file_xprintf(arguments->outfp, "number of special characters: ");
    gt_file_xprintf(arguments->outfp, "%lu (%lu range(s))\n",
                                      gt_encseq_specialcharacters(encseq),
                                      gt_encseq_realspecialranges(encseq));

    gt_file_xprintf(arguments->outfp, "accesstype: ");
    gt_file_xprintf(arguments->outfp, "%s\n",
                   gt_encseq_access_type_str(gt_encseq_accesstype_get(encseq)));

    gt_file_xprintf(arguments->outfp, "bits used per character: ");
    gt_file_xprintf(arguments->outfp, "%f\n",
      (double) ((uint64_t) CHAR_BIT * (uint64_t) gt_encseq_sizeofrep(encseq)) /
      (double) gt_encseq_total_length(encseq));

    gt_file_xprintf(arguments->outfp, "has special ranges: ");
    gt_file_xprintf(arguments->outfp, "%s\n",
                                      gt_encseq_has_specialranges(encseq)
                                        ? "yes"
                                        : "no");

    gt_file_xprintf(arguments->outfp, "has description support: ");
    gt_file_xprintf(arguments->outfp, "%s\n",
                                      gt_encseq_has_description_support(encseq)
                                        ? "yes"
                                        : "no");

    gt_file_xprintf(arguments->outfp, "has multiple sequence support: ");
    gt_file_xprintf(arguments->outfp, "%s\n",
                                      gt_encseq_has_multiseq_support(encseq)
                                        ? "yes"
                                        : "no");
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(encseq_loader);

  return had_err;
}

GtTool* gt_einfo(void)
{
  return gt_tool_new(gt_einfo_arguments_new,
                  gt_einfo_arguments_delete,
                  gt_einfo_option_parser_new,
                  gt_einfo_arguments_check,
                  gt_einfo_runner);
}
