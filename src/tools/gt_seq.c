/*
  Copyright (c) 2006-2010 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/ma.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/undef_api.h"
#include "core/eansi.h"
#include "tools/gt_seq.h"

typedef struct {
  bool recreate,
       showfasta,
       gc_content,
       stat,
       seqlengthdistri,
       seqnum_from_0;
  GtUword showseqnum,
          showseqlength,
                width;
  GtStr *showseqnum_inputfile;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} SeqArguments;

static void* gt_seq_arguments_new(void)
{
  SeqArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  arguments->showseqnum_inputfile = gt_str_new();
  return arguments;
}

static void gt_seq_arguments_delete(void *tool_arguments)
{
  SeqArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_str_delete(arguments->showseqnum_inputfile);
  gt_free(arguments);
}

static GtOptionParser* gt_seq_option_parser_new(void *tool_arguments)
{
  SeqArguments *arguments = tool_arguments;
  GtOption *option, *option_recreate, *option_showfasta, *option_showseqnum,
           *option_width, *option_stat, *option_showseqnum_inputfile,
           *option_seqnum_from_0, *option_showseqlength;
  GtOptionParser *op;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] sequence_file [...]",
                            "Parse the given sequence file(s) and construct "
                            "the corresponding index files.");

  /* -recreate */
  option_recreate = gt_option_new_bool("recreate", "recreate index files, even "
                                       "if they exist already",
                                       &arguments->recreate, false);
  gt_option_parser_add_option(op, option_recreate);

  /* -showfasta */
  option_showfasta = gt_option_new_bool("showfasta", "show all sequences (in "
                                        "FASTA format)", &arguments->showfasta,
                                        false);
  gt_option_parser_add_option(op, option_showfasta);

  /* -showseqnum */
#define GT_COMMENT_COUNTED_FROM_0\
        "(sequences are counted from 1 unless -seqnum_from_0 is specified)"
  option_showseqnum = gt_option_new_uword_min("showseqnum",
                      "show sequence with given number "
                      GT_COMMENT_COUNTED_FROM_0,
                      &arguments->showseqnum,
                      GT_UNDEF_UWORD, 0);
  gt_option_parser_add_option(op, option_showseqnum);

  option_showseqlength = gt_option_new_uword_min("showseqlength",
                      "show length of sequence with given number "
                      GT_COMMENT_COUNTED_FROM_0,
                      &arguments->showseqlength,
                      GT_UNDEF_UWORD, 0);
  gt_option_parser_add_option(op, option_showseqlength);

  /* -showseqnum_inputfile */
  option_showseqnum_inputfile = gt_option_new_filename("showseqnum_inputfile",
                                "show sequence with given numbers "
                                "specified line by line in given file "
                                GT_COMMENT_COUNTED_FROM_0,
                                arguments->showseqnum_inputfile);
  gt_option_parser_add_option(op, option_showseqnum_inputfile);

  /* -seqnum_from_0 */
  option_seqnum_from_0 = gt_option_new_bool("seqnum_from_0",
                        "for options -showseqnum and -showseqnum_inputfile: "
                        "count sequence numbers from 0 instead of 1",
                        &arguments->seqnum_from_0, false);
  gt_option_parser_add_option(op, option_seqnum_from_0);

  /* -gc-content */
  option = gt_option_new_bool("gc-content", "print GC-content (for DNA files)",
                              &arguments->gc_content, false);
  gt_option_parser_add_option(op, option);

  /* -stat */
  option_stat = gt_option_new_bool("stat", "show sequence statistics",
                                   &arguments->stat, false);
  gt_option_parser_add_option(op, option_stat);

  /* -seqlengthdistri */
  option = gt_option_new_bool("seqlengthdistri", "show sequence length "
                              "distribution", &arguments->seqlengthdistri,
                              false);
  gt_option_parser_add_option(op, option);

  /* -width */

  option_width = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, option_width);

  /* output file options */
  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  /* option implications */
  gt_option_imply_either_3(option_width, option_showfasta, option_showseqnum,
                           option_showseqnum_inputfile);
  gt_option_imply_either_3(option_seqnum_from_0, option_showseqnum,
                           option_showseqlength,
                           option_showseqnum_inputfile);

  /* option exclusions */
  gt_option_exclude(option_showfasta, option_stat);
  gt_option_exclude(option_showfasta, option_showseqnum);
  gt_option_exclude(option_showseqnum, option_stat);
  gt_option_exclude(option_showseqnum, option_showseqnum_inputfile);
  gt_option_exclude(option_showseqlength, option_stat);
  gt_option_exclude(option_showseqlength, option_showseqnum);
  gt_option_exclude(option_showseqlength, option_showseqnum_inputfile);

  /* set minimal arguments */
  gt_option_parser_set_min_args(op, 1);

  return op;
}

static int gt_seq_arguments_check(int rest_argc, void *tool_arguments,
                                  GtError *err)
{
  SeqArguments *arguments = tool_arguments;
  gt_error_check(err);
  gt_assert(arguments);
  /* option -showseqnum makes only sense if we got a single sequence file */
  if (arguments->showseqnum != GT_UNDEF_UWORD ||
      arguments->showseqlength != GT_UNDEF_UWORD)
  {
    if (rest_argc > 1)
    {
      gt_error_set(err, "option '-showseq%s' only makes sense with a single "
                        "sequence_file",
                        arguments->showseqnum != GT_UNDEF_UWORD ? "num"
                                                                : "length");
      return -1;
    }
    if ((arguments->showseqnum == 0 || arguments->showseqlength == 0) &&
        !arguments->seqnum_from_0)
    {
      gt_error_set(err, "argument to option '-showseq%s' must be >= 1 if "
                        "option -seqnum_from_0 is not used",
                        arguments->showseqnum == 0 ? "num"
                                                   : "length");
      return -1;
    }
  }
  return 0;
}

static int gt_seq_runner(int argc, const char **argv, int parsed_args,
                         void *tool_arguments, GtError *err)
{
  SeqArguments *arguments = tool_arguments;
  GtBioseq *bioseq;
  int arg = parsed_args, had_err = 0;
  gt_error_check(err);
  gt_assert(tool_arguments);

  while (!had_err && arg < argc) {
    /* bioseq construction */
    if (arguments->recreate)
      bioseq = gt_bioseq_new_recreate(argv[arg], err);
    else
      bioseq = gt_bioseq_new(argv[arg], err);
    if (!bioseq)
      had_err = -1;

    /* output */
    if (!had_err && arguments->showfasta)
      gt_bioseq_show_as_fasta(bioseq, arguments->width, arguments->outfp);

    if (!had_err && (arguments->showseqnum != GT_UNDEF_UWORD ||
                     arguments->showseqlength != GT_UNDEF_UWORD ||
                     gt_str_length(arguments->showseqnum_inputfile) > 0))
    {
      GtUword maxseqnum = gt_bioseq_number_of_sequences(bioseq);
      if (arguments->seqnum_from_0)
      {
        gt_assert(maxseqnum > 0);
        maxseqnum--;
      }
      if (arguments->showseqnum != GT_UNDEF_UWORD ||
          arguments->showseqlength != GT_UNDEF_UWORD)
      {
        if ((arguments->showseqnum != GT_UNDEF_UWORD &&
             arguments->showseqnum > maxseqnum) ||
            (arguments->showseqlength != GT_UNDEF_UWORD &&
             arguments->showseqlength > maxseqnum))
        {
          gt_error_set(err, "argument '" GT_WU "' to option '-showseq%s' "
                            "is too large. The largest possible number is '"
                            GT_WU "'",
                     arguments->showseqnum,
                     arguments->showseqnum != GT_UNDEF_UWORD ? "num"
                                                             : "length",
                     maxseqnum);
          had_err = -1;
        }
        if (!had_err) {
          if (arguments->showseqnum != GT_UNDEF_UWORD)
          {
            gt_bioseq_show_sequence_as_fasta(bioseq,
                                             arguments->seqnum_from_0
                                               ? arguments->showseqnum
                                               : arguments->showseqnum - 1,
                                             arguments->width,
                                             arguments->outfp);
          } else
          {
            GtUword seqnum = arguments->seqnum_from_0
                               ? arguments->showseqlength
                               : arguments->showseqlength - 1;
            gt_file_xprintf(arguments->outfp,GT_WU "\n",
                            gt_bioseq_get_sequence_length(bioseq, seqnum));
          }
        }
      } else
      {
        FILE *fpin;

        gt_assert(gt_str_length(arguments->showseqnum_inputfile) > 0);
        if (strcmp(gt_str_get(arguments->showseqnum_inputfile),"-") == 0)
        {
          fpin = stdin;
        } else
        {
          fpin = gt_efopen(gt_str_get(arguments->showseqnum_inputfile),
                           "rb",err);
        }
        if (fpin == NULL)
        {
          had_err = -1;
        } else
        {
          GtWord seqnum_input;
          while (fscanf(fpin,GT_WD,&seqnum_input) == 1)
          {
            if (seqnum_input < 0 ||
                (!arguments->seqnum_from_0 && seqnum_input == 0) ||
                (GtUword) seqnum_input > maxseqnum)
            {
              gt_error_set(err, "sequence number '" GT_WD "' is not in range "
                                "[%d," GT_WU "]",
                                 seqnum_input,
                                 arguments->seqnum_from_0 ? 0 : 1,
                                 maxseqnum);
              had_err = -1;
            }
            if (!had_err)
            {
              gt_bioseq_show_sequence_as_fasta(bioseq,
                                           arguments->seqnum_from_0
                                             ? (GtUword) seqnum_input
                                             : (GtUword) seqnum_input - 1,
                                           arguments->width, arguments->outfp);
            }
          }
          if (strcmp(gt_str_get(arguments->showseqnum_inputfile),"-") == 0)
          {
            fclose(fpin);
          }
        }
      }
    }

    if (!had_err && arguments->gc_content)
      gt_bioseq_show_gc_content(bioseq, arguments->outfp);

    if (!had_err && arguments->stat)
      gt_bioseq_show_stat(bioseq, arguments->outfp);

    if (!had_err && arguments->seqlengthdistri)
      gt_bioseq_show_seqlengthdistri(bioseq, arguments->outfp);

    /* free */
    gt_bioseq_delete(bioseq);

    arg++;
  }

  return had_err;
}

GtTool* gt_seq(void)
{
  return gt_tool_new(gt_seq_arguments_new,
                     gt_seq_arguments_delete,
                     gt_seq_option_parser_new,
                     gt_seq_arguments_check,
                     gt_seq_runner);
}
