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
#include "core/encseq_metadata.h"
#include "core/ma.h"
#include "core/output_file_api.h"
#include "core/unused_api.h"
#include "tools/gt_encseq_info.h"

typedef struct {
  bool nomap,
       mirror,
       noindexname,
       show_alphabet;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} GtEncseqInfoArguments;

static void* gt_encseq_info_arguments_new(void)
{
  GtEncseqInfoArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_encseq_info_arguments_delete(void *tool_arguments)
{
  GtEncseqInfoArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_encseq_info_option_parser_new(void *tool_arguments)
{
  GtEncseqInfoArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionnomap;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] indexname",
                            "Display meta-information about an "
                            "encoded sequence.");

  optionnomap = gt_option_new_bool("nomap", "do not map encoded sequence "
                                            "(gives less information)",
                                   &arguments->nomap, false);
  gt_option_parser_add_option(op, optionnomap);

  option = gt_option_new_bool("mirrored", "use mirrored encoded sequence "
                                          "(DNA only)",
                              &arguments->mirror, false);
  gt_option_parser_add_option(op, option);
  gt_option_exclude(optionnomap, option);

  option = gt_option_new_bool("noindexname", "do not output index name",
                              &arguments->noindexname, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("show_alphabet", "output alphabet definition",
                              &arguments->show_alphabet, false);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_min_max_args(op, 1, 1);
  return op;
}

static int gt_encseq_info_runner(GT_UNUSED int argc, const char **argv,
                           int parsed_args, void *tool_arguments,
                           GtError *err)
{
  GtEncseqInfoArguments *arguments = tool_arguments;
  int had_err = 0;
  GtAlphabet *alpha;
  const GtUchar *chars;
  gt_error_check(err);
  gt_assert(arguments);

  if (arguments->nomap) {
    GtEncseqMetadata *emd = gt_encseq_metadata_new(argv[parsed_args], err);
    if (!emd)
      had_err = -1;

    if (!had_err) {
      if (!arguments->noindexname) {
        gt_file_xprintf(arguments->outfp, "index name: ");
        gt_file_xprintf(arguments->outfp, "%s\n", argv[parsed_args]);
      }

      gt_file_xprintf(arguments->outfp, "file format version: ");
      gt_file_xprintf(arguments->outfp, "%lu\n",
                                          gt_encseq_metadata_version(emd));

      gt_file_xprintf(arguments->outfp, "64-bit file: ");
      gt_file_xprintf(arguments->outfp, "%s\n", gt_encseq_metadata_is64bit(emd)
                                                  ? "yes"
                                                  : "no");

      gt_file_xprintf(arguments->outfp, "total length: ");
      gt_file_xprintf(arguments->outfp, "%lu\n",
                                        gt_encseq_metadata_total_length(emd));

      gt_file_xprintf(arguments->outfp, "number of sequences: ");
      gt_file_xprintf(arguments->outfp, "%lu\n",
                                      gt_encseq_metadata_num_of_sequences(emd));

      gt_file_xprintf(arguments->outfp, "number of files: ");
      gt_file_xprintf(arguments->outfp, "%lu\n",
                                        gt_encseq_metadata_num_of_files(emd));

      gt_file_xprintf(arguments->outfp, "length of shortest/longest "
                                        "sequence: ");
      gt_file_xprintf(arguments->outfp, "%lu/%lu\n",
                                        gt_encseq_metadata_min_seq_length(emd),
                                        gt_encseq_metadata_max_seq_length(emd));

      gt_file_xprintf(arguments->outfp, "accesstype: ");
      gt_file_xprintf(arguments->outfp, "%s\n",
                 gt_encseq_access_type_str(gt_encseq_metadata_accesstype(emd)));

      alpha = gt_encseq_metadata_alphabet(emd);
      chars = gt_alphabet_characters(alpha);
      gt_file_xprintf(arguments->outfp, "alphabet size: ");
      gt_file_xprintf(arguments->outfp, "%u\n",
                                        gt_alphabet_num_of_chars(alpha));
      gt_file_xprintf(arguments->outfp, "alphabet characters: ");
      gt_file_xprintf(arguments->outfp, "%.*s", gt_alphabet_num_of_chars(alpha),
                                        (char*) chars);
      if (gt_alphabet_is_dna(alpha))
        gt_file_xprintf(arguments->outfp, " (DNA)");
      if (gt_alphabet_is_protein(alpha))
        gt_file_xprintf(arguments->outfp, " (Protein)");
      gt_file_xprintf(arguments->outfp, "\n");
      if (arguments->show_alphabet) {
        GtStr *out = gt_str_new();
        gt_alphabet_to_str(alpha, out);
        gt_file_xprintf(arguments->outfp, "alphabet definition:\n");
        gt_file_xprintf(arguments->outfp, "%s\n", gt_str_get(out));
        gt_str_delete(out);
      }

    }
    gt_encseq_metadata_delete(emd);
  } else {
    GtEncseqLoader *encseq_loader;
    GtEncseq *encseq;

    encseq_loader = gt_encseq_loader_new();
    if (arguments->mirror)
      gt_encseq_loader_mirror(encseq_loader);
    if (!(encseq = gt_encseq_loader_load(encseq_loader,
                                         argv[parsed_args], err)))
      had_err = -1;

    if (!had_err) {
      const GtStrArray *filenames;
      unsigned long i;

      if (!arguments->noindexname) {
        gt_file_xprintf(arguments->outfp, "index name: ");
        gt_file_xprintf(arguments->outfp, "%s\n", argv[parsed_args]);
      }

      gt_file_xprintf(arguments->outfp, "file format version: ");
      gt_file_xprintf(arguments->outfp, "%lu\n", gt_encseq_version(encseq));

      gt_file_xprintf(arguments->outfp, "64-bit file: ");
      gt_file_xprintf(arguments->outfp, "%s\n", gt_encseq_is_64_bit(encseq)
                                                   ? "yes"
                                                   : "no");

      gt_file_xprintf(arguments->outfp, "total length: ");
      gt_file_xprintf(arguments->outfp, "%lu\n",
                                        gt_encseq_total_length(encseq));

      gt_file_xprintf(arguments->outfp, "compressed size: ");
      gt_file_xprintf(arguments->outfp, "%lu bytes\n",
                                        gt_encseq_sizeofrep(encseq));

      gt_file_xprintf(arguments->outfp, "number of sequences: ");
      gt_file_xprintf(arguments->outfp, "%lu\n",
                                        gt_encseq_num_of_sequences(encseq));

      gt_file_xprintf(arguments->outfp, "number of files: ");
      gt_file_xprintf(arguments->outfp, "%lu\n",
                                        gt_encseq_num_of_files(encseq));

      gt_file_xprintf(arguments->outfp, "length of shortest/longest "
                                        "sequence: ");
      gt_file_xprintf(arguments->outfp, "%lu/%lu\n",
                                      gt_encseq_min_seq_length(encseq),
                                      gt_encseq_max_seq_length(encseq));

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
      gt_file_xprintf(arguments->outfp, "%u\n",
                                        gt_alphabet_num_of_chars(alpha));
      gt_file_xprintf(arguments->outfp, "alphabet characters: ");
      gt_file_xprintf(arguments->outfp, "%.*s", gt_alphabet_num_of_chars(alpha),
                                        (char*) chars);
      if (gt_alphabet_is_dna(alpha))
        gt_file_xprintf(arguments->outfp, " (DNA)");
      if (gt_alphabet_is_protein(alpha))
        gt_file_xprintf(arguments->outfp, " (Protein)");
      gt_file_xprintf(arguments->outfp, "\n");
      if (arguments->show_alphabet) {
        GtStr *out = gt_str_new();
        gt_alphabet_to_str(alpha, out);
        gt_file_xprintf(arguments->outfp, "alphabet definition:\n");
        gt_file_xprintf(arguments->outfp, "%s\n", gt_str_get(out));
        gt_str_delete(out);
      }

      gt_file_xprintf(arguments->outfp, "character distribution:\n");
      for (i = 0; i < gt_alphabet_num_of_chars(alpha); i++) {
        unsigned long cc;
        cc = gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, chars[i]));
        gt_file_xprintf(arguments->outfp, "\t%c: %lu (%.2f%%)\n",
                                          (char) chars[i],
                                          cc,
                             (cc /(double) (gt_encseq_total_length(encseq)
                                  - gt_encseq_num_of_sequences(encseq)+1))*100);
      }

      gt_file_xprintf(arguments->outfp, "number of wildcards: ");
      gt_file_xprintf(arguments->outfp, "%lu (%lu range(s))\n",
                                        gt_encseq_wildcards(encseq),
                                        gt_encseq_realwildcardranges(encseq));

      gt_file_xprintf(arguments->outfp, "number of special characters: ");
      gt_file_xprintf(arguments->outfp, "%lu (%lu range(s))\n",
                                        gt_encseq_specialcharacters(encseq),
                                        gt_encseq_realspecialranges(encseq));

      gt_file_xprintf(arguments->outfp, "length of longest non-special "
                                        "character stretch: ");
      gt_file_xprintf(arguments->outfp, "%lu\n",
                                   gt_encseq_lengthoflongestnonspecial(encseq));

      gt_file_xprintf(arguments->outfp, "accesstype: ");
      gt_file_xprintf(arguments->outfp, "%s\n",
                   gt_encseq_access_type_str(gt_encseq_accesstype_get(encseq)));

      gt_file_xprintf(arguments->outfp, "bits used per character: ");
      gt_file_xprintf(arguments->outfp, "%f\n",
        (double) ((uint64_t) CHAR_BIT *
                  (uint64_t) gt_encseq_sizeofrep(encseq)) /
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
  }

  return had_err;
}

GtTool* gt_encseq_info(void)
{
  return gt_tool_new(gt_encseq_info_arguments_new,
                  gt_encseq_info_arguments_delete,
                  gt_encseq_info_option_parser_new,
                  NULL,
                  gt_encseq_info_runner);
}
