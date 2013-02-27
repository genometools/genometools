/*
  Copyright (c) 2008-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2011      Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2008-2011 Center for Bioinformatics, University of Hamburg

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

#include "core/bioseq.h"
#include "core/fasta.h"
#include "core/fileutils_api.h"
#include "core/ma.h"
#include "core/output_file_api.h"
#include "core/option_api.h"
#include "core/progressbar.h"
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/string_distri.h"
#include "core/unused_api.h"
#include "extended/md5set.h"
#include "tools/gt_sequniq.h"

typedef struct {
  bool seqit, verbose, rev;
  unsigned long width, nofseqs;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} GtSequniqArguments;

static void* gt_sequniq_arguments_new(void)
{
  GtSequniqArguments *arguments = gt_calloc((size_t)1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_sequniq_arguments_delete(void *tool_arguments)
{
  GtSequniqArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_sequniq_option_parser_new(void *tool_arguments)
{
  GtSequniqArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *seqit_option, *verbose_option, *width_option, *rev_option,
           *nofseqs_option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] sequence_file [...] ",
                            "Filter out repeated sequences in given in given "
                            "sequence_file(s).");

  /* -seqit */
  seqit_option = gt_option_new_bool("seqit", "use sequence iterator",
                                    &arguments->seqit, false);
  gt_option_is_development_option(seqit_option);
  gt_option_parser_add_option(op, seqit_option);

  /* -nofseqs */
  nofseqs_option = gt_option_new_ulong("nofseqs", "number of sequences "
      "(improves efficiency)\ndefault: unspecified",
      &arguments->nofseqs, 0);
  gt_option_is_development_option(nofseqs_option);
  gt_option_hide_default(nofseqs_option);
  gt_option_parser_add_option(op, nofseqs_option);

  /* -rev */
  rev_option = gt_option_new_bool("rev", "filter out also sequences whose "
      "reverse complement is identical to a sequence already output",
      &arguments->rev, false);
  gt_option_parser_add_option(op, rev_option);

  /* -v */
  verbose_option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, verbose_option);

  /* -width */
  width_option = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, width_option);

  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  /* option implications */
  gt_option_imply(verbose_option, seqit_option);

  gt_option_parser_set_min_args(op, 1U);
  return op;
}

static int gt_sequniq_runner(int argc, const char **argv, int parsed_args,
                             void *tool_arguments, GtError *err)
{
  GtSequniqArguments *arguments = tool_arguments;
  unsigned long long duplicates = 0, num_of_sequences = 0;
  int i, had_err = 0;
  GtMD5Set *md5set;

  gt_error_check(err);
  gt_assert(arguments);
  md5set = gt_md5set_new(arguments->nofseqs);
  if (!arguments->seqit) {
    unsigned long j;
    GtBioseq *bs;

    for (i = parsed_args; !had_err && i < argc; i++) {
      if (!(bs = gt_bioseq_new(argv[i], err)))
        had_err = -1;
      if (!had_err) {
        GtMD5SetStatus retval;
        for (j = 0; j < gt_bioseq_number_of_sequences(bs) && !had_err; j++) {
          char *seq = gt_bioseq_get_sequence(bs, j);
          retval = gt_md5set_add_sequence(md5set, seq,
                                          gt_bioseq_get_sequence_length(bs, j),
                                          arguments->rev, err);
          if (retval == GT_MD5SET_NOT_FOUND)
            gt_fasta_show_entry(gt_bioseq_get_description(bs, j), seq,
                                gt_bioseq_get_sequence_length(bs, j),
                                arguments->width, arguments->outfp);
          else if (retval != GT_MD5SET_ERROR)
            duplicates++;
          else
            had_err = -1;
          num_of_sequences++;
          gt_free(seq);
        }
        gt_bioseq_delete(bs);
      }
    }
  }
  else {
    GtSeqIterator *seqit;
    GtStrArray *files;
    off_t totalsize;
    const GtUchar *sequence;
    char *desc;
    unsigned long len;

    files = gt_str_array_new();
    for (i = parsed_args; i < argc; i++)
      gt_str_array_add_cstr(files, argv[i]);
    totalsize = gt_files_estimate_total_size(files);
    seqit = gt_seq_iterator_sequence_buffer_new(files, err);
    if (!seqit)
      had_err = -1;
    if (!had_err) {
      if (arguments->verbose) {
        gt_progressbar_start(gt_seq_iterator_getcurrentcounter(seqit,
                                                            (unsigned long long)
                                                            totalsize),
                             (unsigned long long) totalsize);
      }
      while (!had_err) {
        GtMD5SetStatus retval;
        if ((gt_seq_iterator_next(seqit, &sequence, &len, &desc, err)) != 1)
          break;

        retval = gt_md5set_add_sequence(md5set, (const char*) sequence, len,
            arguments->rev, err);
        if (retval == GT_MD5SET_NOT_FOUND)
          gt_fasta_show_entry(desc, (const char*) sequence, len,
                              arguments->width, arguments->outfp);
        else if (retval != GT_MD5SET_ERROR)
          duplicates++;
        else
          had_err = -1;
        num_of_sequences++;
      }
      if (arguments->verbose)
        gt_progressbar_stop();
      gt_seq_iterator_delete(seqit);
    }
    gt_str_array_delete(files);
  }

  /* show statistics */
  if (!had_err) {
    fprintf(stderr, "# %lu out of %lu sequences have been removed (%.3f%%)\n",
            (unsigned long)duplicates, (unsigned long)num_of_sequences,
            ((double) duplicates / (double)num_of_sequences) * 100.0);
  }

  gt_md5set_delete(md5set);
  return had_err;
}

GtTool* gt_sequniq(void)
{
  return gt_tool_new(gt_sequniq_arguments_new,
                     gt_sequniq_arguments_delete,
                     gt_sequniq_option_parser_new,
                     NULL,
                     gt_sequniq_runner);
}
