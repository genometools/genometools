/*
  Copyright (c) 2009-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2009-2013 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies
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

#include "core/fileutils_api.h"
#include "core/filelengthvalues.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/sequence_buffer.h"
#include "core/seq_iterator_sequence_buffer.h"
#include "core/versionfunc.h"
#include "core/progressbar.h"
#include "extended/reverse_api.h"
#include "core/unused_api.h"
#include "tools/gt_convertseq.h"

typedef struct
{
  bool reduce_wc_dna,
       reduce_wc_prot,
       revcomp,
       showflv,
       showseq,
       verbose;
  GtUword fastawidth;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} GtConvertseqArguments;

static void* gt_convertseq_arguments_new(void)
{
  GtConvertseqArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  arguments->outfp = NULL;
  return arguments;
}

static void gt_convertseq_arguments_delete(void *tool_arguments)
{
  GtConvertseqArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_convertseq_option_parser_new(void *tool_arguments)
{
  GtConvertseqArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *o;
  gt_assert(arguments);

  op = gt_option_parser_new("[options] file [...]",
                            "Parse and convert sequence file formats "
                            "(FASTA/FASTQ, GenBank, EMBL).");

  o = gt_option_new_bool("v","be verbose", &arguments->verbose, false);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_bool("r","reverse complement sequences",
                         &arguments->revcomp, false);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_bool("showfilelengthvalues","show filelengths",
                         &arguments->showflv, false);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_bool("noseq","do not show sequences",
                         &arguments->showseq, false);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_uword("fastawidth",
                          "FASTA output line width, 0 for unlimited",
                         &arguments->fastawidth, 60UL);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_bool("contractdnawc", "replace stretches of DNA wildcards "
                                          "with a single 'N'",
                         &arguments->reduce_wc_dna, false);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_bool("contractproteinwc", "replace stretches of protein "
                                              "wildcards with a single 'X'",
                         &arguments->reduce_wc_prot, false);
  gt_option_parser_add_option(op, o);

  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  gt_option_parser_set_min_args(op, 1U);

  return op;
}

static int gt_convertseq_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, GtError *err)
{
  GtConvertseqArguments *arguments = tool_arguments;
  int had_err = 0, i;
  GtFilelengthvalues *flv;
  GtSeqIterator *seqit;
  GtSequenceBuffer *sb = NULL;
  GtStrArray *files;
  const GtUchar *sequence;
  char *desc;
  GtUword len, j;
  off_t totalsize;
  gt_error_check(err);
  gt_assert(arguments != NULL);

  files = gt_str_array_new();
  for (i = parsed_args; i < argc; i++)
  {
    gt_str_array_add_cstr(files, argv[i]);
  }
  totalsize = gt_files_estimate_total_size(files);

  flv = gt_calloc((size_t) gt_str_array_size(files),
                  sizeof (GtFilelengthvalues));

  sb = gt_sequence_buffer_new_guess_type(files, err);
  if (!sb) {
    had_err = -1;
  }
  if (!had_err) {
    gt_sequence_buffer_set_filelengthtab(sb, flv);
    /* read input using seqiterator */
    seqit = gt_seq_iterator_sequence_buffer_new_with_buffer(sb);
    if (arguments->verbose)
    {
      gt_progressbar_start(gt_seq_iterator_getcurrentcounter(seqit,
                                                           (GtUint64)
                                                           totalsize),
                           (GtUint64) totalsize);
    }
    while (true)
    {
      GtUchar *seq = NULL;
      desc = NULL;
      j = 0UL;
      had_err = gt_seq_iterator_next(seqit, &sequence, &len, &desc, err);
      if (had_err != 1)
        break;
      if (arguments->revcomp) {
        GtUchar *newseq = gt_calloc((size_t) len+1, sizeof (GtUchar));
        memcpy(newseq, sequence, (size_t) len*sizeof (GtUchar));
        had_err = gt_reverse_complement((char*) newseq, len, err);
        if (had_err)
          break;
        seq = newseq;
      } else seq = (GtUchar*) sequence;

      if (!arguments->showseq) {
        bool in_wildcard = false;
        gt_file_xprintf(arguments->outfp, ">%s\n", desc);
        for (i = 0; (GtUword) i < len; i++) {
          if (arguments->reduce_wc_dna) {
            switch (seq[i]) {
              case 'a':
              case 'A':
              case 'c':
              case 'C':
              case 'g':
              case 'G':
              case 't':
              case 'u':
              case 'T':
              case 'U':
                in_wildcard = false;
                gt_file_xfputc((int) seq[i], arguments->outfp);
                j++;
                break;
              default:
                if (!in_wildcard) {
                  in_wildcard = true;
                  if (isupper((int) seq[i]))
                    gt_file_xfputc((int) 'N', arguments->outfp);
                  else
                    gt_file_xfputc((int) 'n', arguments->outfp);
                  j++;
                }
            }
          }
          else if (arguments->reduce_wc_prot) {
            switch (seq[i]) {
              case 'X':
              case 'B':
              case 'Z':
                if (!in_wildcard) {
                  in_wildcard = true;
                  gt_file_xfputc((int) 'N', arguments->outfp);
                  j++;
                }
                break;
              case 'x':
              case 'b':
              case 'z':
                if (!in_wildcard) {
                  in_wildcard = true;
                  gt_file_xfputc((int) 'n', arguments->outfp);
                  j++;
                }
                break;
              default:
                in_wildcard = false;
                gt_file_xfputc((int) seq[i], arguments->outfp);
                j++;
            }
          }
          else {
            gt_file_xfputc((int) seq[i], arguments->outfp);
            j++;
          }
          if (arguments->fastawidth > 0 && j % arguments->fastawidth == 0) {
            j = 0;
            gt_file_xprintf(arguments->outfp, "\n");
          }
        }
        if (arguments->fastawidth == 0 || len % arguments->fastawidth != 0)
            gt_file_xprintf(arguments->outfp, "\n");
      }
      if (arguments->revcomp) {
        gt_free(seq);
      }
    }
    if (arguments->showflv) {
      for (j=0;j<gt_str_array_size(files);j++) {
        fprintf(stderr, "file "GT_WU" (%s): "GT_WU"/"GT_WU"\n",
               j,
               gt_str_array_get(files, j),
               (GtUword) flv[j].length,
               (GtUword) flv[j].effectivelength);
      }
    }
    if (arguments->verbose)
    {
      gt_progressbar_stop();
    }
    gt_sequence_buffer_delete(sb);
    gt_seq_iterator_delete(seqit);
  }
  gt_str_array_delete(files);
  gt_free(flv);

  return had_err;
}

GtTool* gt_convertseq(void)
{
  return gt_tool_new(gt_convertseq_arguments_new,
                     gt_convertseq_arguments_delete,
                     gt_convertseq_option_parser_new,
                     NULL,
                     gt_convertseq_runner);
}
