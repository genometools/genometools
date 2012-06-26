/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/undef_api.h"
#include "core/seqiterator_fastq.h"
#include "core/unused_api.h"
#include "extended/feature_type.h"
#include "extended/hpol_processor.h"
#include "tools/gt_hopcorrect.h"

typedef struct {
  GtStr  *encseqinput, *annotation, *map, *outfilename;
  unsigned long hmin, max_hlen_diff;
  bool verbose, map_is_bam, rchk, stats;
  double minc;
  GtStrArray *readset;
} GtHopcorrectArguments;

static void* gt_hopcorrect_arguments_new(void)
{
  GtHopcorrectArguments *arguments = gt_malloc(sizeof (*arguments));
  arguments->encseqinput = gt_str_new();
  arguments->annotation = gt_str_new();
  arguments->map = gt_str_new();
  arguments->outfilename = gt_str_new();
  arguments->readset = gt_str_array_new();
  return arguments;
}

static void gt_hopcorrect_arguments_delete(void *tool_arguments)
{
  GtHopcorrectArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->encseqinput);
  gt_str_delete(arguments->annotation);
  gt_str_delete(arguments->map);
  gt_str_delete(arguments->outfilename);
  gt_str_array_delete(arguments->readset);
  gt_free(arguments);
}

static GtOptionParser* gt_hopcorrect_option_parser_new(void *tool_arguments)
{
  GtHopcorrectArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("-r <encseq> [-m <sam|bam>] "
      "[-a <gff3>] [-o <outfile>]",
      "Reference-based homopolymer length adjustment.");

  /* -r */
  option = gt_option_new_string("r", "reference sequence\n"
      "format: GtEncseq",
      arguments->encseqinput, NULL);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* -a */
  option = gt_option_new_string("a", "annotation of reference sequence\n"
      "it must be sorted by coordinates on the reference sequence\n"
      "(this can be e.g. done using: gt gff3 -sort)\n"
      "format: sorted GFF3",
      arguments->annotation, NULL);
  gt_option_parser_add_option(op, option);

  /* -m */
  option = gt_option_new_string("m", "mapping over reference sequence\n"
      "it must be sorted by coordinates on the reference sequence\n"
      "(this can be e.g. done using: samtools sort)\n"
      "format: sorted SAM/BAM",
      arguments->map, NULL);
  gt_option_parser_add_option(op, option);

  /* -rchk */
  option = gt_option_new_bool("rchk", "check that ref region of aligned\n"
      "segments is identical to encseq",
      &arguments->rchk, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -o */
  option = gt_option_new_string("o", "output file for corrected sequences\n"
      "format: FASTQ",
      arguments->outfilename, NULL);
 gt_option_parser_add_option(op, option);

  /* -bam */
  option = gt_option_new_bool("bam", "set to true if map is BAM, "
      "false if SAM\ndefault: true (BAM)",
      &arguments->map_is_bam, true);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* -hmin */
  option = gt_option_new_ulong_min("hmin", "minimal homopolymer length",
      &arguments->hmin, 3UL, 2UL);
  gt_option_parser_add_option(op, option);

  /* -maxd */
  option = gt_option_new_ulong("maxd", "maximal difference "
      "in homopolymer length by which a correction is applied\n"
      "default: infinite", &arguments->max_hlen_diff, GT_UNDEF_ULONG);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* -minc */
#ifndef S_SPLINT_S
  option = gt_option_new_double_min_max("minc", "mimimal homopolymer length "
      "consensus in segments to a different value than the length in the "
      "reference, after which no correction to the segments is applied;\n"
      "expressed as decimal value between 0.0 and 1.0;\ne.g. 0.9 means: do not "
      "correct any segment if the length of the homopolymer in 90\% or more "
      "of the segments agrees on a different value than the length in the "
      "reference\ndefault: undefined", &arguments->minc, 2.0D, 0, 1.0D);
#endif
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* -stats */
  option = gt_option_new_bool("stats", "output statistics for each "
      "correction position", &arguments->stats, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -outorder */
  option = gt_option_new_filename_array("outorder",
      "specify the name of the files (FastQ format) containing the "
      "uncorrected readset; this allows to output the reads in "
      "the same order, increasing memory usage", arguments->readset);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_max_args(op, 0, 0);

  return op;
}

static int gt_hopcorrect_runner(GT_UNUSED int argc, GT_UNUSED const char **argv,
    GT_UNUSED int parsed_args, void *tool_arguments, GtError *err)
{
  GtHopcorrectArguments *arguments = tool_arguments;
  GtEncseq *encseq;
  GtEncseqLoader *el;
  GtLogger *v_logger;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);
  v_logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  gt_assert(gt_str_length(arguments->encseqinput) > 0);
  el = gt_encseq_loader_new();
  gt_encseq_loader_drop_description_support(el);
  gt_encseq_loader_disable_autosupport(el);
  encseq = gt_encseq_loader_load(el, gt_str_get(arguments->encseqinput), err);
  if (encseq == NULL)
  {
    had_err = -1;
  }
  if (!had_err)
  {
    GtHpolProcessor *hpp;
    GtSeqposClassifier *spc = NULL;
    GtSamfileIterator *sfi = NULL;
    GtAlignedSegmentsPile *asp = NULL;
    GtFile *outfile = NULL;
    GtAlphabet *dna = gt_alphabet_new_dna();
    GtSeqIterator *readset_iter = NULL;
    hpp = gt_hpol_processor_new(encseq, arguments->hmin);
    if (gt_str_length(arguments->map) > 0)
    {
      if (arguments->map_is_bam)
        sfi = gt_samfile_iterator_new_bam(gt_str_get(arguments->map), dna, err);
      else
        sfi = gt_samfile_iterator_new_sam(gt_str_get(arguments->map), dna,
            NULL, err);
      if (sfi == NULL)
        had_err = -1;
      else
      {
        asp = gt_aligned_segments_pile_new(sfi);
        if (!arguments->rchk)
          gt_hpol_processor_enable_segments_hlen_adjustment(hpp, asp,
              arguments->max_hlen_diff, arguments->minc);
        else
          gt_hpol_processor_enable_aligned_segments_refregionscheck(hpp, asp);
        if (arguments->stats)
          gt_hpol_processor_enable_statistics_output(hpp, NULL);
        if (!had_err && gt_str_length(arguments->outfilename) > 0)
        {
          outfile = gt_file_new(gt_str_get(arguments->outfilename), "w", err);
          if (outfile == NULL)
            had_err = -1;
          else
            gt_hpol_processor_enable_segments_output(hpp, outfile);
          if (!had_err && gt_str_array_size(arguments->readset) > 0)
          {
            readset_iter = gt_seqiterator_fastq_new(arguments->readset, err);
            if (readset_iter == NULL)
              had_err = -1;
            else
              gt_hpol_processor_sort_segments_output(hpp, readset_iter);
          }
        }
      }
    }
    if (!had_err && gt_str_length(arguments->annotation) > 0)
    {
      spc = gt_seqpos_classifier_new(gt_str_get(arguments->annotation),
          gt_ft_CDS);
      gt_hpol_processor_restrict_to_feature_type(hpp, spc);
    }
    if (!had_err)
      had_err = gt_hpol_processor_run(hpp, v_logger, err);
    gt_aligned_segments_pile_delete(asp);
    gt_samfile_iterator_delete(sfi);
    gt_seqpos_classifier_delete(spc);
    gt_file_delete(outfile);
    gt_hpol_processor_delete(hpp);
    gt_seqiterator_delete(readset_iter);
    gt_alphabet_delete(dna);
  }
  gt_logger_delete(v_logger);
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  return had_err;
}

GtTool* gt_hopcorrect(void)
{
  return gt_tool_new(gt_hopcorrect_arguments_new,
                  gt_hopcorrect_arguments_delete,
                  gt_hopcorrect_option_parser_new,
                  NULL,
                  gt_hopcorrect_runner);
}
