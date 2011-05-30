/*
  Copyright (c) 2008-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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
#include "core/fa.h"
#include "core/fasta.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/outputfile.h"
#include "core/string_distri.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/gtdatahelp.h"
#include "tools/gt_fingerprint.h"

typedef struct {
  bool show_duplicates;
  GtStr *checklist,
        *extract;
  unsigned long width;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} FingerprintArguments;

static void* gt_fingerprint_arguments_new(void)
{
  FingerprintArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->checklist = gt_str_new();
  arguments->extract = gt_str_new();
  arguments->ofi = gt_outputfileinfo_new();
  return arguments;
}

static void gt_fingerprint_arguments_delete(void *tool_arguments)
{
  FingerprintArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_outputfileinfo_delete(arguments->ofi);
  gt_str_delete(arguments->extract);
  gt_str_delete(arguments->checklist);
  gt_free(arguments);
}

static GtOptionParser* gt_fingerprint_option_parser_new(GT_UNUSED
                                                        void *tool_arguments)
{
  FingerprintArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *check_option, *duplicates_option, *extract_option, *width_option;
  gt_assert(arguments);
  op = gt_option_parser_new("[option ...] sequence_file [...] ",
                            "Compute MD5 fingerprints for each sequence given "
                            "in sequence_file(s).");

  /* -check */
  check_option = gt_option_new_filename("check", "compare all fingerprints "
                                        "contained in the given checklist file "
                                        "with checksums in given "
                                        "sequence_files(s). The comparison is "
                                        "successful, if all fingerprints given "
                                        "in checkfile can be found in the "
                                        "sequence_file(s) in the exact same "
                                        "quantity and vice versa.",
                                        arguments->checklist);
  gt_option_parser_add_option(op, check_option);

  /* -duplicates */
  duplicates_option = gt_option_new_bool("duplicates", "show duplicate "
                                         "fingerprints from given "
                                         "sequence_file(s).",
                                         &arguments->show_duplicates, false);
  gt_option_parser_add_option(op, duplicates_option);

  /* -extract */
  extract_option = gt_option_new_string("extract",
                                        "extract the sequence(s) with "
                                        "the given fingerprint from "
                                        "sequence_file(s) and show them on "
                                        "stdout.", arguments->extract, NULL);
  gt_option_parser_add_option(op, extract_option);

  /* -width */
  width_option = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, width_option);

  gt_outputfile_register_options(op, &arguments->outfp, arguments->ofi);

  /* option exclusions */
  gt_option_exclude(check_option, duplicates_option);
  gt_option_exclude(extract_option, check_option);
  gt_option_exclude(extract_option, duplicates_option);

  /* option implications */
  gt_option_imply(width_option, extract_option);

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  gt_option_parser_set_min_args(op, 1);
  return op;
}

static void proc_superfluous_sequence(const char *string,
                                      GT_UNUSED unsigned long occurrences,
                                      GT_UNUSED double probability, void *data)
{
  bool *comparisons_failed = data;
  gt_assert(string && occurrences && comparisons_failed);
  printf("%s only in sequence_file(s)\n", string);
  *comparisons_failed = true;
}

static int compare_fingerprints(GtStringDistri *sd, const char *checklist,
                                GtError *err)
{
  bool comparisons_failed = false, use_stdin = false;
  FILE *checkfile;
  GtStr *line;
  gt_error_check(err);
  gt_assert(sd && checklist);
  if (!strcmp(checklist, "-"))
    use_stdin = true;
  checkfile = use_stdin ? stdin : gt_fa_xfopen(checklist, "r");
  line = gt_str_new();
  /* process checklist */
  while (gt_str_read_next_line(line, checkfile) != EOF) {
    if (gt_string_distri_get(sd, gt_str_get(line)))
      gt_string_distri_sub(sd, gt_str_get(line));
    else {
      printf("%s only in checklist\n", gt_str_get(line));
      comparisons_failed = true;
    }
    gt_str_reset(line);
  }
  gt_str_delete(line);
  if (!use_stdin)
    gt_fa_xfclose(checkfile);
  /* process remaining sequence_file(s) fingerprints */
  gt_string_distri_foreach(sd, proc_superfluous_sequence, &comparisons_failed);
  if (comparisons_failed) {
    gt_error_set(err, "fingerprint comparison failed");
    return -1;
  }
  return 0;
}

typedef struct {
  unsigned long long duplicates,
                     num_of_sequences;
} FingerprintInfo;

static void show_duplicate(const char *fingerprint, unsigned long occurrences,
                           GT_UNUSED double probability, void *data)
{
  FingerprintInfo *info = data;
  if (occurrences > 1) {
    printf("%s\t%lu\n", fingerprint, occurrences);
    info->duplicates += occurrences - 1;
  }
  info->num_of_sequences += occurrences;
}

static int show_duplicates(GtStringDistri *sd, GtError *err)
{
  FingerprintInfo info;
  gt_error_check(err);
  gt_assert(sd);
  info.duplicates = 0;
  info.num_of_sequences = 0;
  gt_string_distri_foreach(sd, show_duplicate, &info);
  if (info.duplicates) {
    gt_error_set(err, "duplicates found: %llu out of %llu (%.3f%%)",
                 info.duplicates, info.num_of_sequences,
                 (((double) info.duplicates / info.num_of_sequences) * 100.0));
    return -1;
  }
  return 0;
}

static int gt_fingerprint_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  FingerprintArguments *arguments = tool_arguments;
  GtBioseq *bs;
  GtStringDistri *sd;
  unsigned long i, j;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);
  sd = gt_string_distri_new();

  /* process sequence files */
  for (i = parsed_args; !had_err && i < argc; i++) {
    if (!(bs = gt_bioseq_new(argv[i], err)))
      had_err = -1;
    if (!had_err) {
      for (j = 0; j < gt_bioseq_number_of_sequences(bs); j++) {
        if (gt_str_length(arguments->checklist) || arguments->show_duplicates)
          gt_string_distri_add(sd, gt_bioseq_get_md5_fingerprint(bs, j));
        else if (gt_str_length(arguments->extract)) {
          if (!strcmp(gt_bioseq_get_md5_fingerprint(bs, j),
                      gt_str_get(arguments->extract))) {
            gt_fasta_show_entry(gt_bioseq_get_description(bs, j),
                                gt_bioseq_get_sequence(bs, j),
                                gt_bioseq_get_sequence_length(bs, j),
                                arguments->width, arguments->outfp);
          }
        }
        else
          gt_xputs(gt_bioseq_get_md5_fingerprint(bs, j));
      }
    }
    gt_bioseq_delete(bs);
  }

  if (!had_err) {
    if (gt_str_length(arguments->checklist))
      had_err = compare_fingerprints(sd, gt_str_get(arguments->checklist), err);
    else if (arguments->show_duplicates)
      had_err = show_duplicates(sd, err);
  }

  gt_string_distri_delete(sd);

  return had_err;
}

GtTool* gt_fingerprint(void)
{
  return gt_tool_new(gt_fingerprint_arguments_new,
                     gt_fingerprint_arguments_delete,
                     gt_fingerprint_option_parser_new,
                     NULL,
                     gt_fingerprint_runner);
}
