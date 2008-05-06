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
#include "libgtcore/fa.h"
#include "libgtcore/fasta.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/string_distri.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"
#include "libgtext/gtdatahelp.h"
#include "tools/gt_fingerprint.h"

typedef struct {
  bool show_duplicates;
  Str *checklist,
      *extract;
} FingerprintArguments;

static void* gt_fingerprint_arguments_new(void)
{
  FingerprintArguments *arguments = ma_calloc(1, sizeof *arguments);
  arguments->checklist = str_new();
  arguments->extract = str_new();
  return arguments;
}

static void gt_fingerprint_arguments_delete(void *tool_arguments)
{
  FingerprintArguments *arguments = tool_arguments;
  if (!arguments) return;
  str_delete(arguments->extract);
  str_delete(arguments->checklist);
  ma_free(arguments);
}

static OptionParser* gt_fingerprint_option_parser_new(UNUSED
                                                      void *tool_arguments)
{
  FingerprintArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *check_option, *duplicates_option, *extract_option;
  assert(arguments);
  op = option_parser_new("[option ...] sequence_file [...] ",
                         "Compute MD5 fingerprints for each sequence given "
                         "in sequence_file(s).");

  /* -check */
  check_option = option_new_filename("check", "Compare all fingerprints "
                                     "contained in the given checklist file "
                                     "with checksums in given "
                                     "sequence_files(s). The comparison is "
                                     "successful, if all fingerprints given in "
                                     "checkfile can be found in the "
                                     "sequence_file(s) in the exact same "
                                     "quantity and vice versa.",
                                     arguments->checklist);
  option_parser_add_option(op, check_option);

  /* -duplicates */
  duplicates_option = option_new_bool("duplicates", "Show duplicate "
                                      "fingerprints from given "
                                      "sequence_file(s).",
                                      &arguments->show_duplicates, false);
  option_parser_add_option(op, duplicates_option);

  /* -extract */
  extract_option = option_new_string("extract", "Extract the sequence(s) with "
                                     "the given fingerprint from "
                                     "sequence_file(s) and show them on "
                                     "stdout.", arguments->extract, NULL);
  option_parser_add_option(op, extract_option);

  /* option exclusions */
  option_exclude(check_option, duplicates_option);
  option_exclude(extract_option, check_option);
  option_exclude(extract_option, duplicates_option);

  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  option_parser_set_min_args(op, 1);
  return op;
}

static void proc_superfluous_sequence(const char *string,
                                      UNUSED unsigned long occurrences,
                                      UNUSED double probability, void *data)
{
  bool *comparisons_failed = data;
  assert(string && occurrences && comparisons_failed);
  printf("%s only in sequence_file(s)\n", string);
  *comparisons_failed = true;
}

static int compare_fingerprints(StringDistri *sd, const char *checklist,
                                Error *err)
{
  bool comparisons_failed = false, use_stdin = false;
  FILE *checkfile;
  Str *line;
  error_check(err);
  assert(sd && checklist);
  if (!strcmp(checklist, "-"))
    use_stdin = true;
  checkfile = use_stdin ? stdin : fa_xfopen(checklist, "r");
  line = str_new();
  /* process checklist */
  while (str_read_next_line(line, checkfile) != EOF) {
    if (string_distri_get(sd, str_get(line)))
      string_distri_sub(sd, str_get(line));
    else {
      printf("%s only in checklist\n", str_get(line));
      comparisons_failed = true;
    }
    str_reset(line);
  }
  str_delete(line);
  if (!use_stdin)
    fa_xfclose(checkfile);
  /* process remaining sequence_file(s) fingerprints */
  string_distri_foreach(sd, proc_superfluous_sequence, &comparisons_failed);
  if (comparisons_failed) {
    error_set(err, "fingerprint comparison failed");
    return -1;
  }
  return 0;
}

typedef struct {
  unsigned long long duplicates,
                     num_of_sequences;
} FingerprintInfo;

static void show_duplicate(const char *fingerprint, unsigned long occurrences,
                           UNUSED double probability, void *data)
{
  FingerprintInfo *info = data;
  if (occurrences > 1) {
    printf("%s\t%lu\n", fingerprint, occurrences);
    info->duplicates += occurrences - 1;
  }
  info->num_of_sequences += occurrences;
}

static int show_duplicates(StringDistri *sd, Error *err)
{
  FingerprintInfo info;
  error_check(err);
  assert(sd);
  info.duplicates = 0;
  info.num_of_sequences = 0;
  string_distri_foreach(sd, show_duplicate, &info);
  if (info.duplicates) {
    error_set(err, "duplicates found: %llu out of %llu (%.3f%%)",
              info.duplicates, info.num_of_sequences,
              (((double) info.duplicates / info.num_of_sequences) * 100.0));
    return -1;
  }
  return 0;
}

static int gt_fingerprint_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, Error *err)
{
  FingerprintArguments *arguments = tool_arguments;
  Bioseq *bs;
  StringDistri *sd;
  unsigned long i, j;
  int had_err = 0;

  error_check(err);
  assert(arguments);
  sd = string_distri_new();

  /* process sequence files */
  for (i = parsed_args; !had_err && i < argc; i++) {
    if (!(bs = bioseq_new(argv[i], err)))
      had_err = -1;
    if (!had_err) {
      for (j = 0; j < bioseq_number_of_sequences(bs); j++) {
        if (str_length(arguments->checklist) || arguments->show_duplicates)
          string_distri_add(sd, bioseq_get_md5_fingerprint(bs, j));
        else if (str_length(arguments->extract)) {
          if (!strcmp(bioseq_get_md5_fingerprint(bs, j),
                      str_get(arguments->extract))) {
            fasta_show_entry(bioseq_get_description(bs, j),
                             bioseq_get_sequence(bs, j),
                             bioseq_get_sequence_length(bs, j), 0);
          }
        }
        else
          xputs(bioseq_get_md5_fingerprint(bs, j));
      }
    }
    bioseq_delete(bs);
  }

  if (!had_err) {
    if (str_length(arguments->checklist))
      had_err = compare_fingerprints(sd, str_get(arguments->checklist), err);
    else if (arguments->show_duplicates)
      had_err = show_duplicates(sd, err);
  }

  string_distri_delete(sd);

  return had_err;
}

Tool* gt_fingerprint(void)
{
  return tool_new(gt_fingerprint_arguments_new,
                  gt_fingerprint_arguments_delete,
                  gt_fingerprint_option_parser_new,
                  NULL,
                  gt_fingerprint_runner);
}
