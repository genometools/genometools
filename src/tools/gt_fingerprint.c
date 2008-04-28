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
#include "libgtcore/stringdistri.h"
#include "libgtcore/unused.h"
#include "libgtext/gtdatahelp.h"
#include "tools/gt_fingerprint.h"

typedef struct {
  bool check,
       show_duplicates;
  Str *checklist;
} FingerprintArguments;

static void* gt_fingerprint_arguments_new(void)
{
  FingerprintArguments *arguments = ma_calloc(1, sizeof *arguments);
  arguments->checklist = str_new();
  return arguments;
}

static void gt_fingerprint_arguments_delete(void *tool_arguments)
{
  FingerprintArguments *arguments = tool_arguments;
  if (!arguments) return;
  str_delete(arguments->checklist);
  ma_free(arguments);
}

static OptionParser* gt_fingerprint_option_parser_new(UNUSED
                                                      void *tool_arguments)
{
  FingerprintArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *check_option, *duplicates_option;
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

  /* option exclusions */
  option_exclude(check_option, duplicates_option);

  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  option_parser_set_min_args(op, 1);
  return op;
}

typedef struct {
  unsigned long long duplicates,
                     num_of_sequences;
} FingerprintInfo;

static void show_duplicates(const char *fingerprint, unsigned long occurrences,
                            UNUSED double probability, void *data)
{
  FingerprintInfo *info = data;
  if (occurrences > 1) {
    printf("%s %lu\n", fingerprint, occurrences);
    info->duplicates += occurrences - 1;
  }
  info->num_of_sequences += occurrences;
}

static int gt_fingerprint_runner(int argc, const char **argv,
                                 void *tool_arguments, Error *err)
{
  FingerprintArguments *arguments = tool_arguments;
  Bioseq *bs;
  StringDistri *sd;
  unsigned long i, j;
  FingerprintInfo info;
  int had_err = 0;

  error_check(err);
  assert(arguments);
  sd = stringdistri_new();

  for (i = 0; !had_err && i < argc; i++) {
    if (!(bs = bioseq_new(argv[i], err)))
      had_err = -1;
    if (!had_err) {
      for (j = 0; j < bioseq_number_of_sequences(bs); j++)
        stringdistri_add(sd, bioseq_get_md5_fingerprint(bs, j));
    }
    bioseq_delete(bs);
  }

  info.duplicates = 0;
  info.num_of_sequences = 0;

  if (!had_err && arguments->show_duplicates) {
    stringdistri_foreach(sd, show_duplicates, &info);
    if (info.duplicates) {
      printf("total number of duplicates: %llu out of %llu (%.3f%%)\n",
             info.duplicates, info.num_of_sequences,
             (((double) info.duplicates / info.num_of_sequences) * 100.0));
    }
  }

  stringdistri_delete(sd);

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
