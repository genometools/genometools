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
#include "libgtcore/fasta.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/ma.h"
#include "libgtcore/md5_fingerprint.h"
#include "libgtcore/outputfile.h"
#include "libgtcore/option.h"
#include "libgtcore/progressbar.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/string_distri.h"
#include "libgtcore/unused.h"
#include "tools/gt_sequniq.h"

typedef struct {
  bool seqit,
       verbose;
  OutputFileInfo *ofi;
  GenFile *outfp;
} SequniqArguments;

static void* gt_sequniq_arguments_new(void)
{
  SequniqArguments *arguments = ma_calloc(1, sizeof *arguments);
  arguments->ofi = outputfileinfo_new();
  return arguments;
}

static void gt_sequniq_arguments_delete(void *tool_arguments)
{
  SequniqArguments *arguments = tool_arguments;
  if (!arguments) return;
  genfile_close(arguments->outfp);
  outputfileinfo_delete(arguments->ofi);
  ma_free(arguments);
}

static OptionParser* gt_sequniq_option_parser_new(void *tool_arguments)
{
  SequniqArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *seqit_option, *verbose_option;
  assert(arguments);

  op = option_parser_new("[option ...] sequence_file [...] ",
                         "Filter out repeated sequences in given in given "
                         "sequence_file(s).");

  /* -seqit */
  seqit_option = option_new_bool("seqit", "use sequence iterator",
                                 &arguments->seqit, false);
  option_is_development_option(seqit_option);
  option_parser_add_option(op, seqit_option);

  /* -v */
  verbose_option = option_new_verbose(&arguments->verbose);
  option_is_development_option(verbose_option);
  option_parser_add_option(op, verbose_option);

  /* option implications */
  option_imply(verbose_option, seqit_option);

  outputfile_register_options(op, &arguments->outfp, arguments->ofi);
  option_parser_set_min_args(op, 1);
  return op;
}

static int gt_sequniq_runner(int argc, const char **argv, int parsed_args,
                             void *tool_arguments, Error *err)
{
  SequniqArguments *arguments = tool_arguments;
  Bioseq *bs;
  StringDistri *sd;
  unsigned long long duplicates = 0, num_of_sequences = 0;
  StrArray *files;
  int had_err = 0;
  SeqIterator *seqit;
  const Uchar *sequence;
  char *desc;
  unsigned long len;
  off_t totalsize;

  error_check(err);
  assert(arguments);
  sd = string_distri_new();

  if (!arguments->seqit) {
    unsigned long i, j;
    for (i = parsed_args; !had_err && i < argc; i++) {
      if (!(bs = bioseq_new(argv[i], err)))
        had_err = -1;
      if (!had_err) {
        for (j = 0; j < bioseq_number_of_sequences(bs); j++) {
          if (!string_distri_get(sd, bioseq_get_md5_fingerprint(bs, j))) {
            string_distri_add(sd, bioseq_get_md5_fingerprint(bs, j));
            fasta_show_entry_generic(bioseq_get_description(bs, j),
                                     bioseq_get_sequence(bs, j),
                                     bioseq_get_sequence_length(bs, j), 0,
                                     arguments->outfp);
          }
          else
            duplicates++;
          num_of_sequences++;
        }
      }
      bioseq_delete(bs);
    }
  }
  else {
    int i;
    files = strarray_new();
    for (i = parsed_args; i < argc; i++)
      strarray_add_cstr(files, argv[i]);
    totalsize = files_estimate_total_size(files);
    seqit = seqiterator_new(files, NULL, true);
    if (arguments->verbose) {
      progressbar_start(seqiterator_getcurrentcounter(seqit,
                                                      (unsigned long long)
                                                      totalsize),
                                                      (unsigned long long)
                                                      totalsize);
    }
    for (;;) {
      char *md5;
      if ((seqiterator_next(seqit, &sequence, &len, &desc, err)) != 1)
        break;
      md5 = md5_fingerprint((const char*) sequence, (unsigned long) len);
      if (!string_distri_get(sd, md5)) {
        string_distri_add(sd, md5);
        fasta_show_entry_generic(desc, (const char*) sequence, len, 0,
                                 arguments->outfp);
      }
      else
        duplicates++;
      num_of_sequences++;
      ma_free(desc);
      ma_free(md5);
    }
    if (arguments->verbose)
      progressbar_stop();
    seqiterator_delete(seqit);
    strarray_delete(files);
  }

  /* show statistics */
  if (!had_err) {
    fprintf(stderr, "# %llu out of %llu sequences have been removed (%.3f%%)\n",
            duplicates, num_of_sequences,
            ((double) duplicates / num_of_sequences) * 100.0);
  }
  string_distri_delete(sd);

  return had_err;
}

Tool* gt_sequniq(void)
{
  return tool_new(gt_sequniq_arguments_new,
                  gt_sequniq_arguments_delete,
                  gt_sequniq_option_parser_new,
                  NULL,
                  gt_sequniq_runner);
}
