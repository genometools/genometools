/*
  Copyright (c) 2008-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <ctype.h>
#include <string.h>
#include "core/bioseq.h"
#include "core/bioseq_col.h"
#include "core/fa.h"
#include "core/fasta.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/seq_info_cache.h"
#include "core/string_distri.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/gtdatahelp.h"
#include "tools/gt_fingerprint.h"

typedef struct {
  bool show_duplicates,
       detect_collisions;
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
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_fingerprint_arguments_delete(void *tool_arguments)
{
  FingerprintArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_str_delete(arguments->extract);
  gt_str_delete(arguments->checklist);
  gt_free(arguments);
}

static GtOptionParser* gt_fingerprint_option_parser_new(GT_UNUSED
                                                        void *tool_arguments)
{
  FingerprintArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *check_option, *collisions_option, *duplicates_option,
           *extract_option, *width_option;
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

  /* -collisions */
  collisions_option = gt_option_new_bool("collisions", "detect hash collisions",
                                         &arguments->detect_collisions, false);
  gt_option_is_development_option(collisions_option);
  gt_option_parser_add_option(op, collisions_option);

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

  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  /* option exclusions */
  gt_option_exclude(check_option, duplicates_option);
  gt_option_exclude(extract_option, check_option);
  gt_option_exclude(extract_option, duplicates_option);
  gt_option_exclude(collisions_option, check_option);
  gt_option_exclude(collisions_option, duplicates_option);
  gt_option_exclude(collisions_option, extract_option);

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

static int compare_md5s(GtBioseqCol *bsc, const GtSeqInfo *si,
                        unsigned long filenum, unsigned long seqnum,
                        const char *md5, GtError *err)
{
  unsigned long i, seq_a_len, seq_b_len;
  char *seq_a_upper, *seq_b_upper;
  char *seq_a, *seq_b;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(bsc && si && md5);
  seq_a_len = gt_seq_col_get_sequence_length((GtSeqCol*) bsc, si->filenum,
                                             si->seqnum);
  seq_b_len = gt_seq_col_get_sequence_length((GtSeqCol*) bsc, filenum, seqnum);
  seq_a = gt_seq_col_get_sequence((GtSeqCol*) bsc, si->filenum, si->seqnum, 0,
                                  seq_a_len - 1);
  seq_b = gt_seq_col_get_sequence((GtSeqCol*) bsc, filenum, seqnum, 0,
                                  seq_b_len - 1);
  seq_a_upper = gt_malloc((seq_a_len + 1) * sizeof (char));
  seq_b_upper = gt_malloc((seq_b_len + 1) * sizeof (char));
  for (i = 0; i < seq_a_len; i++)
    seq_a_upper[i] = toupper(seq_a[i]);
  for (i = 0; i < seq_b_len; i++)
    seq_b_upper[i] = toupper(seq_b[i]);
  gt_free(seq_a);
  gt_free(seq_b);
  seq_a_upper[seq_a_len] = '\0';
  seq_b_upper[seq_b_len] = '\0';
  if (strcmp(seq_a_upper, seq_b_upper)) {
    gt_error_set(err, "sequence collision detected for fingerprint '%s'", md5);
    had_err = -1;
  }
  gt_free(seq_b_upper);
  gt_free(seq_a_upper);
  return had_err;
}

static int detect_collisions_on_bsc(GtBioseqCol *bsc, GtError *err)
{
  unsigned long filenum, seqnum;
  GtSeqInfoCache *sic;
  gt_error_check(err);
  gt_assert(bsc);
  int had_err = 0;
  sic = gt_seq_info_cache_new();
  for (filenum = 0;
       !had_err && filenum < gt_seq_col_num_of_files((GtSeqCol*) bsc);
       filenum++) {
    for (seqnum = 0;
         !had_err && seqnum < gt_seq_col_num_of_seqs((GtSeqCol*) bsc, filenum);
         seqnum++) {
      const GtSeqInfo *si_ptr;
      const char *md5;
      md5 = gt_seq_col_get_md5_fingerprint((GtSeqCol*) bsc, filenum, seqnum);
      if ((si_ptr = gt_seq_info_cache_get(sic, md5)))
        had_err = compare_md5s(bsc, si_ptr, filenum, seqnum, md5, err);
      else {
        GtSeqInfo si;
        si.filenum = filenum;
        si.seqnum = seqnum;
        gt_seq_info_cache_add(sic, md5, &si);
      }
    }
  }
  gt_seq_info_cache_delete(sic);
  return had_err;
}

static int detect_collisions(int num_of_seqfiles, const char **seqfiles,
                             GtError *err)
{
  GtStrArray *sequence_files;
  GtSeqCol *sc;
  unsigned long i;
  int had_err = 0;
  gt_error_check(err);
  sequence_files = gt_str_array_new();
  for (i = 0; i < num_of_seqfiles; i++)
    gt_str_array_add_cstr(sequence_files, seqfiles[i]);
  if (!(sc = gt_bioseq_col_new(sequence_files, err)))
    had_err = -1;
  if (!had_err)
    had_err = detect_collisions_on_bsc((GtBioseqCol*) sc, err);
  gt_seq_col_delete(sc);
  gt_str_array_delete(sequence_files);
  return had_err;
}

static int gt_fingerprint_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  FingerprintArguments *arguments = tool_arguments;
  bool extract_found = true;
  GtBioseq *bs;
  GtStringDistri *sd;
  unsigned long i, j;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);
  sd = gt_string_distri_new();

  if (gt_str_length(arguments->extract))
    extract_found = false;

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
            char *seq = gt_bioseq_get_sequence(bs, j);
            gt_fasta_show_entry(gt_bioseq_get_description(bs, j),
                                seq,
                                gt_bioseq_get_sequence_length(bs, j),
                                arguments->width, arguments->outfp);
            gt_free(seq);
            extract_found = true;
          }
        }
        else if (!arguments->detect_collisions)
          gt_xputs(gt_bioseq_get_md5_fingerprint(bs, j));
      }
    }
    gt_bioseq_delete(bs);
  }

  if (!had_err && !extract_found) {
    gt_error_set(err, "could not find sequence with fingerprint '%s' in given "
                      "sequence file(s)", gt_str_get(arguments->extract));
    had_err = -1;
  }

  if (!had_err) {
    if (gt_str_length(arguments->checklist))
      had_err = compare_fingerprints(sd, gt_str_get(arguments->checklist), err);
    else if (arguments->show_duplicates)
      had_err = show_duplicates(sd, err);
    else if (arguments->detect_collisions)
      had_err = detect_collisions(argc - parsed_args, argv + parsed_args, err);
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
