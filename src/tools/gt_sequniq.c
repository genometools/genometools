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
#include <ctype.h>
#include "md5.h"
#include "core/bioseq.h"
#include "core/fasta.h"
#include "core/fileutils_api.h"
#include "core/hashtable.h"
#include "core/ma.h"
#include "core/outputfile.h"
#include "core/option_api.h"
#include "core/progressbar.h"
#include "core/safearith.h"
#include "core/seqiterator_sequence_buffer.h"
#include "core/string_distri.h"
#include "core/unused_api.h"
#include "extended/reverse.h"
#include "tools/gt_sequniq.h"

typedef struct {
  bool seqit,
       verbose,
       r;
  unsigned long width;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
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
  gt_file_delete(arguments->outfp);
  gt_outputfileinfo_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_sequniq_option_parser_new(void *tool_arguments)
{
  GtSequniqArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *seqit_option, *verbose_option, *width_option, *r_option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] sequence_file [...] ",
                            "Filter out repeated sequences in given in given "
                            "sequence_file(s).");

  /* -seqit */
  seqit_option = gt_option_new_bool("seqit", "use sequence iterator",
                                    &arguments->seqit, false);
  gt_option_is_development_option(seqit_option);
  gt_option_parser_add_option(op, seqit_option);

  /* -r */
  r_option = gt_option_new_bool("r", "filter out also sequences whose"
      "reverse complement is identical to sequence already output",
      &arguments->r, false);
  gt_option_is_development_option(r_option);
  gt_option_imply(r_option, seqit_option);
  gt_option_parser_add_option(op, r_option);

  /* -v */
  verbose_option = gt_option_new_verbose(&arguments->verbose);
  gt_option_is_development_option(verbose_option);
  gt_option_parser_add_option(op, verbose_option);

  /* -width */
  width_option = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, width_option);

  gt_outputfile_register_options(op, &arguments->outfp, arguments->ofi);

  /* option implications */
  gt_option_imply(verbose_option, seqit_option);

  gt_option_parser_set_min_args(op, 1);
  return op;
}

htsize_t gt_md5set_md5_hash(const void *md5)
{
  return gt_uint32_data_hash(md5, sizeof (uint32_t) / sizeof (unsigned char));
}

int gt_md5set_md5_cmp(const void *md5A, const void *md5B)
{
  return memcmp(md5A, md5B, sizeof (unsigned char) * 16);
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
  const GtUchar *sequence;
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
            gt_fasta_show_entry(gt_bioseq_get_description(bs, j),
                                gt_bioseq_get_sequence(bs, j),
                                gt_bioseq_get_sequence_length(bs, j),
                                arguments->width, arguments->outfp);
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
    char *revcompl = NULL, *upper = NULL;
    unsigned long revcompl_alloc = 0;
    GtHashtable *md5set;
    static const HashElemInfo hashtype =
        {gt_md5set_md5_hash, {NULL}, 16, gt_md5set_md5_cmp, NULL, NULL};
    char md5hash[16], md5hash_rc[16];
    unsigned long j;

    md5set = gt_hashtable_new(hashtype);
    files = gt_str_array_new();
    for (i = parsed_args; i < argc; i++)
      gt_str_array_add_cstr(files, argv[i]);
    totalsize = gt_files_estimate_total_size(files);
    seqit = gt_seqiterator_sequence_buffer_new(files, err);
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
        if ((gt_seqiterator_next(seqit, &sequence, &len, &desc, err)) != 1)
          break;

        /* check size of necessary buffers */
        if (upper == NULL)
          upper = gt_malloc(sizeof (char) * len);
        else if (revcompl_alloc < len)
          upper = gt_realloc(upper, sizeof (char) * len);
        if (arguments->r)
        {
          if (revcompl == NULL)
            revcompl = gt_malloc(sizeof (char) * len);
          else if (revcompl_alloc < len + 1)
            revcompl = gt_realloc(revcompl, sizeof (char) * len);
        }

        for (j = 0; j < len; j++)
          upper[j] = toupper(sequence[j]);

        md5(upper, gt_safe_cast2long(len), md5hash);

        if (arguments->r)
        {
          (void)memcpy(revcompl, upper, len);
          had_err = gt_reverse_complement(revcompl, len, err);
          if (had_err)
            break;
          md5(revcompl, gt_safe_cast2long(len), md5hash_rc);
        }
        if (!gt_hashtable_get(md5set, md5hash) &&
            (!arguments->r || !gt_hashtable_get(md5set, md5hash_rc)))
        {
          gt_hashtable_add(md5set, md5hash);
          if (arguments->r)
            gt_hashtable_add(md5set, md5hash_rc);
          gt_fasta_show_entry(desc, (const char*) sequence, len,
                              arguments->width, arguments->outfp);
        }
        else
          duplicates++;
        num_of_sequences++;
      }
      gt_free(upper);
      if (arguments->r)
        gt_free(revcompl);
      if (arguments->verbose)
        gt_progressbar_stop();
      gt_seqiterator_delete(seqit);
    }
    gt_hashtable_delete(md5set);
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
