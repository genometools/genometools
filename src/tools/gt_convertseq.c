/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2009-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef S_SPLINT_S
#include <ctype.h>
#else
#ifndef GT_ISUPPER_DEFINED
#define GT_ISUPPER_DEFINED
int isupper(int c);
#endif
#endif
#include <string.h>
#include "core/fileutils_api.h"
#include "core/filelengthvalues.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/sequence_buffer.h"
#include "core/seq_iterator_sequence_buffer.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "core/progressbar.h"
#include "extended/reverse_api.h"
#include "tools/gt_convertseq.h"

typedef struct
{
  bool reduce_wc_dna,
       reduce_wc_prot,
       revcomp,
       showflv,
       showseq,
       verbose;
  unsigned long fastawidth;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} ConvertseqOptions;

static GtOPrval parse_options(ConvertseqOptions *opts, int *parsed_args,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *optionverbose, *o;
  GtOPrval oprval;

  gt_error_check(err);

  op = gt_option_parser_new("[options] file [...]",
                            "Parse and convert sequence file formats.");

  optionverbose = gt_option_new_bool("v","be verbose",
                                     &opts->verbose, false);
  gt_option_parser_add_option(op, optionverbose);

  optionverbose = gt_option_new_bool("r","reverse complement sequences",
                                     &opts->revcomp, false);
  gt_option_parser_add_option(op, optionverbose);

  o = gt_option_new_bool("showfilelengthvalues","show filelengths",
                         &opts->showflv, false);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_bool("noseq","do not show sequences",
                         &opts->showseq, false);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_ulong("fastawidth","FASTA output line width",
                         &opts->fastawidth, 60UL);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_bool("contractdnawc", "replace stretches of DNA wildcards "
                                          "with a single 'N'",
                         &opts->reduce_wc_dna, false);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_bool("contractproteinwc", "replace stretches of protein "
                                              "wildcards with a single 'X'",
                         &opts->reduce_wc_prot, false);
  gt_option_parser_add_option(op, o);

  gt_output_file_register_options(op, &opts->outfp, opts->ofi);

  gt_option_parser_set_min_args(op, 1U);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_convertseq(int argc, const char **argv, GtError *err)
{
  GtStrArray *files;
  const GtUchar *sequence;
  char *desc;
  unsigned long len, j;
  int i, parsed_args, had_err = 0;
  off_t totalsize;
  ConvertseqOptions opts;
  GtFilelengthvalues *flv;
  GtSeqIterator *seqit;
  GtSequenceBuffer *sb = NULL;
  opts.ofi = gt_output_file_info_new();
  opts.outfp = NULL;

  gt_error_check(err);

  /* option parsing */
  switch (parse_options(&opts,&parsed_args, argc, argv, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR:
        gt_file_delete(opts.outfp);
        gt_output_file_info_delete(opts.ofi);
        return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT:
        gt_file_delete(opts.outfp);
        gt_output_file_info_delete(opts.ofi);
        return 0;
  }

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
    if (opts.verbose)
    {
      gt_progressbar_start(gt_seq_iterator_getcurrentcounter(seqit,
                                                           (unsigned long long)
                                                           totalsize),
                           (unsigned long long) totalsize);
    }
    while (true)
    {
      GtUchar *seq = NULL;
      desc = NULL;
      i = 0; j = 0UL;
      had_err = gt_seq_iterator_next(seqit, &sequence, &len, &desc, err);
      if (had_err != 1)
        break;
      if (opts.revcomp) {
        GtUchar *newseq = gt_calloc((size_t) len+1, sizeof (GtUchar));
        memcpy(newseq, sequence, (size_t) len*sizeof (GtUchar));
        had_err = gt_reverse_complement((char*) newseq, len, err);
        if (had_err)
          break;
        seq = newseq;
      } else seq = (GtUchar*) sequence;

      if (!opts.showseq) {
        bool in_wildcard = false;
        gt_file_xprintf(opts.outfp, ">%s\n", desc);
        for (i = 0; (unsigned long) i < len; i++) {
          if (opts.reduce_wc_dna) {
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
                gt_file_xfputc((int) seq[i], opts.outfp);
                j++;
                break;
              default:
                if (!in_wildcard) {
                  in_wildcard = true;
                  if (isupper((int) seq[i]))
                    gt_file_xfputc((int) 'N', opts.outfp);
                  else
                    gt_file_xfputc((int) 'n', opts.outfp);
                  j++;
                }
            }
          }
          else if (opts.reduce_wc_prot) {
            switch (seq[i]) {
              case 'X':
              case 'B':
              case 'Z':
                if (!in_wildcard) {
                  in_wildcard = true;
                  gt_file_xfputc((int) 'N', opts.outfp);
                  j++;
                }
                break;
              case 'x':
              case 'b':
              case 'z':
                if (!in_wildcard) {
                  in_wildcard = true;
                  gt_file_xfputc((int) 'n', opts.outfp);
                  j++;
                }
                break;
              default:
                in_wildcard = false;
                gt_file_xfputc((int) seq[i], opts.outfp);
                j++;
            }
          }
          else {
            gt_file_xfputc((int) seq[i], opts.outfp);
            j++;
          }
          if ((j % opts.fastawidth) == 0) {
            j = 0;
            gt_file_xprintf(opts.outfp, "\n");
          }
        }
        gt_file_xprintf(opts.outfp, "\n");
      }
      if (opts.revcomp) {
        gt_free(seq);
      }
    }
    if (opts.showflv) {
      for (j=0;j<gt_str_array_size(files);j++) {
        fprintf(stderr, "file %lu (%s): %lu/%lu\n",
               j,
               gt_str_array_get(files, j),
               (unsigned long) flv[j].length,
               (unsigned long) flv[j].effectivelength);
      }
    }
    gt_sequence_buffer_delete(sb);
    gt_seq_iterator_delete(seqit);
    if (opts.verbose)
    {
      gt_progressbar_stop();
    }
  }
  gt_file_delete(opts.outfp);
  gt_output_file_info_delete(opts.ofi);
  gt_str_array_delete(files);
  gt_free(flv);
  return had_err;
}
