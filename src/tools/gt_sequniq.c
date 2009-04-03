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
#include "core/fasta.h"
#include "core/fileutils.h"
#include "core/ma.h"
#include "core/md5_fingerprint.h"
#include "core/outputfile.h"
#include "core/option.h"
#include "core/progressbar.h"
#include "core/seqiterator.h"
#include "core/string_distri.h"
#include "core/unused_api.h"
#include "tools/gt_sequniq.h"

typedef struct {
  bool seqit,
       verbose;
  GtOutputFileInfo *ofi;
  GtGenFile *outfp;
} GtSequniqArguments;

static void* gt_sequniq_arguments_new(void)
{
  GtSequniqArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_outputfileinfo_new();
  return arguments;
}

static void gt_sequniq_arguments_delete(void *tool_arguments)
{
  GtSequniqArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_genfile_close(arguments->outfp);
  gt_outputfileinfo_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_sequniq_option_parser_new(void *tool_arguments)
{
  GtSequniqArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *seqit_option, *verbose_option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] sequence_file [...] ",
                         "Filter out repeated sequences in given in given "
                         "sequence_file(s).");

  /* -seqit */
  seqit_option = gt_option_new_bool("seqit", "use sequence iterator",
                                 &arguments->seqit, false);
  gt_option_is_development_option(seqit_option);
  gt_option_parser_add_option(op, seqit_option);

  /* -v */
  verbose_option = gt_option_new_verbose(&arguments->verbose);
  gt_option_is_development_option(verbose_option);
  gt_option_parser_add_option(op, verbose_option);

  /* option implications */
  gt_option_imply(verbose_option, seqit_option);

  gt_outputfile_register_options(op, &arguments->outfp, arguments->ofi);
  gt_option_parser_set_min_args(op, 1);
  return op;
}

static int gt_sequniq_runner(int argc, const char **argv, int parsed_args,
                             void *tool_arguments, GtError *err)
{
  GtSequniqArguments *arguments = tool_arguments;
  GtBioseq *bs;
  GtStringDistri *sd;
  unsigned long long duplicates = 0, num_of_sequences = 0;
  GtStrArray *files;
  int had_err = 0;
  GtSeqIterator *seqit;
  const Uchar *sequence;
  char *desc;
  unsigned long len;
  off_t totalsize;

  gt_error_check(err);
  gt_assert(arguments);
  sd = gt_string_distri_new();

  if (!arguments->seqit) {
    unsigned long i, j;
    for (i = parsed_args; !had_err && i < argc; i++) {
      if (!(bs = gt_bioseq_new(argv[i], err)))
        had_err = -1;
      if (!had_err) {
        for (j = 0; j < gt_bioseq_number_of_sequences(bs); j++) {
          if (!gt_string_distri_get(sd, gt_bioseq_get_md5_fingerprint(bs, j))) {
            gt_string_distri_add(sd, gt_bioseq_get_md5_fingerprint(bs, j));
            gt_fasta_show_entry_generic(gt_bioseq_get_description(bs, j),
                                        gt_bioseq_get_sequence(bs, j),
                                        gt_bioseq_get_sequence_length(bs, j), 0,
                                        arguments->outfp);
          }
          else
            duplicates++;
          num_of_sequences++;
        }
      }
      gt_bioseq_delete(bs);
    }
  }
  else {
    int i;
    files = gt_str_array_new();
    for (i = parsed_args; i < argc; i++)
      gt_str_array_add_cstr(files, argv[i]);
    totalsize = gt_files_estimate_total_size(files);
    seqit = gt_seqiterator_new(files, err);
    if (!seqit)
      had_err = -1;
    if (!had_err) {
      if (arguments->verbose) {
        gt_progressbar_start(gt_seqiterator_getcurrentcounter(seqit,
                                                            (unsigned long long)
                                                            totalsize),
                             (unsigned long long) totalsize);
      }
      for (;;) {
        char *md5;
        if ((gt_seqiterator_next(seqit, &sequence, &len, &desc, err)) != 1)
          break;
        md5 = gt_md5_fingerprint((const char*) sequence, (unsigned long) len);
        if (!gt_string_distri_get(sd, md5)) {
          gt_string_distri_add(sd, md5);
          gt_fasta_show_entry_generic(desc, (const char*) sequence, len, 0,
                                      arguments->outfp);
        }
        else
          duplicates++;
        num_of_sequences++;
        gt_free(desc);
        gt_free(md5);
      }
      if (arguments->verbose)
        gt_progressbar_stop();
      gt_seqiterator_delete(seqit);
    }
    gt_str_array_delete(files);
  }

  /* show statistics */
  if (!had_err) {
    fprintf(stderr, "# %llu out of %llu sequences have been removed (%.3f%%)\n",
            duplicates, num_of_sequences,
            ((double) duplicates / num_of_sequences) * 100.0);
  }
  gt_string_distri_delete(sd);

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
