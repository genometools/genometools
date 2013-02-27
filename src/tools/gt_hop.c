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

#include "core/basename_api.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "core/seq_iterator_fastq_api.h"
#include "core/unused_api.h"
#include "extended/feature_type.h"
#include "extended/gtdatahelp.h"
#include "extended/hpol_processor.h"
#include "tools/gt_hop.h"

typedef struct {
  GtStr  *encseqinput, *annotation, *map, *outfilename, *atype, *outprefix;
  unsigned long hmin, clenmax, read_hmin, covmin, mapqmin;
  bool verbose, map_is_sam, rchk, stats, allow_partial, allow_multiple,
       aggressive, moderate, conservative, expert, state_of_truth;
  double altmax, refmin;
  GtStrArray *readset;
} GtHopArguments;

static void* gt_hop_arguments_new(void)
{
  GtHopArguments *arguments = gt_malloc(sizeof (*arguments));
  arguments->encseqinput = gt_str_new();
  arguments->annotation = gt_str_new();
  arguments->atype = gt_str_new();
  arguments->map = gt_str_new();
  arguments->outfilename = gt_str_new();
  arguments->outprefix = gt_str_new();
  arguments->readset = gt_str_array_new();
  return arguments;
}

static void gt_hop_arguments_delete(void *tool_arguments)
{
  GtHopArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->encseqinput);
  gt_str_delete(arguments->annotation);
  gt_str_delete(arguments->atype);
  gt_str_delete(arguments->map);
  gt_str_delete(arguments->outfilename);
  gt_str_delete(arguments->outprefix);
  gt_str_array_delete(arguments->readset);
  gt_free(arguments);
}

static GtOptionParser* gt_hop_option_parser_new(void *tool_arguments)
{
  GtHopArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *ann_option, *aggressive_option, *conservative_option,
           *moderate_option, *expert_option, *o_option, *reads_option,
           *stats_option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("-<mode> -ref <encseq> -map <sam/bam> "
      "-reads <fastq> [options...]",
      "Reference-based homopolymer error correction.");

  /* -ref */
  option = gt_option_new_string("ref",
      "reference sequence in GtEncseq format\n"
      "(can be prepared using gt encseq encode)",
      arguments->encseqinput, NULL);
  gt_option_is_mandatory(option);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* -map */
  option = gt_option_new_string("map",
      "mapping of reads to reference\n"
      "it must be in SAM/BAM format, and sorted by coordinate\n"
      "(can be prepared e.g. using: samtools sort)",
      arguments->map, NULL);
  gt_option_is_mandatory(option);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* -sam */
  option = gt_option_new_bool("sam",
      "mapping file is SAM\ndefault: BAM",
      &arguments->map_is_sam, false);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* -aggressive */
  aggressive_option = gt_option_new_bool("aggressive",
      "correct as much as possible",
       &arguments->aggressive, false);
  gt_option_hide_default(aggressive_option);
  gt_option_parser_add_option(op, aggressive_option);

  /* -moderate */
  moderate_option = gt_option_new_bool("moderate",
      "mediate between sensitivity and precision",
      &arguments->moderate, false);
  gt_option_hide_default(moderate_option);
  gt_option_exclude(moderate_option, aggressive_option);
  gt_option_parser_add_option(op, moderate_option);

  /* -conservative */
  conservative_option = gt_option_new_bool("conservative",
      "correct only most likely errors",
      &arguments->conservative, false);
  gt_option_hide_default(conservative_option);
  gt_option_exclude(conservative_option, aggressive_option);
  gt_option_exclude(conservative_option, moderate_option);
  gt_option_parser_add_option(op, conservative_option);

  /* -expert */
  expert_option = gt_option_new_bool("expert",
      "manually select correction criteria",
      &arguments->expert, false);
  gt_option_hide_default(expert_option);
  gt_option_exclude(expert_option, aggressive_option);
  gt_option_exclude(expert_option, moderate_option);
  gt_option_exclude(expert_option, conservative_option);
  gt_option_parser_add_option(op, expert_option);

  /* -reads */
  reads_option = gt_option_new_filename_array("reads",
      "uncorrected read file(s) in FastQ format;\n"
      "the corrected reads are output in the currect working directory "
      "in files which are named as the input files, each prepended "
      "by a prefix (see -outprefix option)\n"
      "-reads allows one to output the reads in the same order as in the input "
      "and is mandatory if the SAM contains more than a single primary "
      "alignment for each read (e.g. output of bwasw)\n"
      "see also -o option as an alternative",
      arguments->readset);
  gt_option_parser_add_option(op, reads_option);

  /* -outprefix */
  option = gt_option_new_string("outprefix",
      "prefix for output filenames (corrected reads)"
      "when -reads is specified\n"
      "the prefix is prepended to each input filename",
      arguments->outprefix, "hop_");
  gt_option_imply(option, reads_option);
  gt_option_parser_add_option(op, option);

  /* -o */
  o_option = gt_option_new_string("o", "output file for corrected reads\n"
      "(see also -reads/-outprefix) if -o is used, reads are output "
      "in a single file in the order they are found in the SAM file "
      "(which usually differ from the original order)\n"
      "this will only work if the reads were aligned with a software which "
      "only includes 1 alignment for each read (e.g. bwa)",
      arguments->outfilename, NULL);
  gt_option_exclude(reads_option, o_option);
  gt_option_is_mandatory_either(reads_option, o_option);
  gt_option_parser_add_option(op, o_option);

  /* -hmin */
  option = gt_option_new_ulong_min("hmin",
      "minimal homopolymer length in reference\n"
      "minimal number of consecutive identical symbols on the reference",
      &arguments->hmin, 3UL, 2UL);
  gt_option_is_extended_option(option);
  gt_option_imply(option, expert_option);
  gt_option_parser_add_option(op, option);

  /* -read-hmin */
  option = gt_option_new_ulong_min("read-hmin",
      "minimal homopolymer length in reads",
      &arguments->read_hmin, 2UL, 1UL);
  gt_option_is_extended_option(option);
  gt_option_imply(option, expert_option);
  gt_option_parser_add_option(op, option);

  /* -altmax */
  option = gt_option_new_double_min_max("altmax",
      "max support of alternate homopol. length;\n"
      "e.g. 0.8 means: do not correct any read if homop. length in more than "
      "80\% of the reads has the same value, different from the reference\n"
      "if altmax is set to 1.0 reads are always corrected",
      &arguments->altmax, (double) 0.8, 0.0, (double) 1.0);
  gt_option_is_extended_option(option);
  gt_option_imply(option, expert_option);
  gt_option_parser_add_option(op, option);

  /* -refmin */
  option = gt_option_new_double_min_max("refmin",
      "min support of reference homopol. length;\n"
      "e.g. 0.1 means: do not correct any read if ref. homop. length "
      "is not present in at least 10\% of the reads\n"
      "if refmin is set to 0.0 reads are always corrected",
      &arguments->refmin, (double) 0.1, 0.0, (double) 1.0);
  gt_option_hide_default(option);
  gt_option_is_extended_option(option);
  gt_option_imply(option, expert_option);
  gt_option_parser_add_option(op, option);

  /* -mapqmin */
  option = gt_option_new_ulong("mapqmin",
      "minimal mapping quality",
      &arguments->mapqmin, 21UL);
  gt_option_is_extended_option(option);
  gt_option_imply(option, expert_option);
  gt_option_parser_add_option(op, option);

  /* -covmin */
  option = gt_option_new_ulong_min("covmin",
      "minimal coverage;\n"
      "e.g. 5 means: do not correct any read if coverage "
      "(number of reads mapped over whole homopolymer) "
      "is less than 5\n"
      "if covmin is set to 1 reads are always corrected",
      &arguments->covmin, 1UL, 1UL);
  gt_option_is_extended_option(option);
  gt_option_imply(option, expert_option);
  gt_option_parser_add_option(op, option);

  /* -allow-multiple */
  option = gt_option_new_bool("allow-muliple",
      "allow multiple corrections in a read",
      &arguments->allow_multiple, false);
  gt_option_is_extended_option(option);
  gt_option_imply(option, expert_option);
  gt_option_parser_add_option(op, option);

  /* -clenmax */
  option = gt_option_new_ulong("clenmax", "maximal correction length\n"
      "default: unlimited", &arguments->clenmax, GT_UNDEF_ULONG);
  gt_option_hide_default(option);
  gt_option_is_extended_option(option);
  gt_option_imply(option, expert_option);
  gt_option_parser_add_option(op, option);

  /* -ann */
  ann_option = gt_option_new_string("ann",
      "annotation of reference sequence\n"
      "it must be sorted by coordinates on the reference sequence\n"
      "(this can be e.g. done using: gt gff3 -sort)\n"
      "if -ann is used, corrections will be limited to homopolymers starting"
      "or ending inside the feature type indicated by -ft option"
      "format: sorted GFF3",
      arguments->annotation, NULL);
  gt_option_is_extended_option(ann_option);
  gt_option_parser_add_option(op, ann_option);

  /* -ft */
  option = gt_option_new_string("ft", "feature type to use when -ann option "
      "is specified",
      arguments->atype, gt_ft_CDS);
  gt_option_is_extended_option(option);
  gt_option_imply(option, ann_option);
  gt_option_parser_add_option(op, option);

  /* -stats */
  stats_option = gt_option_new_bool("stats", "output statistics for each "
      "correction position", &arguments->stats, false);
  gt_option_is_development_option(stats_option);
  gt_option_parser_add_option(op, stats_option);

  /* -state-of-truth */
  option = gt_option_new_bool("state-of-truth",
      "similar to -stats in -aggressive "
      "mode, but used to determine the \"state of truth\" set of corrections "
      "for evaluation; currently the only difference is that if multiple hits "
      "are present for a read, they are used all independently for "
      "correction (-reads must be set)",
      &arguments->state_of_truth, false);
  gt_option_is_development_option(option);
  gt_option_exclude(option, stats_option);
  gt_option_exclude(option, aggressive_option);
  gt_option_exclude(option, moderate_option);
  gt_option_exclude(option, conservative_option);
  gt_option_exclude(option, expert_option);
  gt_option_parser_add_option(op, option);

  /* -rchk */
  option = gt_option_new_bool("rchk", "debug option; check that ref region "
      "of aligned segments is compatible with encseq data",
      &arguments->rchk, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -allow-partial */
  option = gt_option_new_bool("allow-partial",
      "allow insertions also if there are less gaps in read homopolymer "
      "than the difference in length with the reference\n"
      "(at most as many symbols as the gaps will be inserted)",
      &arguments->allow_partial, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);

  gt_option_parser_set_min_max_args(op, 0, 0);

  return op;
}

int gt_hop_arguments_check(GT_UNUSED int rest_argc, void *tool_arguments,
    GtError* err)
{
  GtHopArguments *args = (GtHopArguments*) tool_arguments;
  int had_err = 0;

  if (!args->aggressive &&
      !args->moderate &&
      !args->conservative &&
      !args->expert &&
      !args->state_of_truth)
  {
    gt_error_set(err, "Select correction mode: "
        "-aggressive, -moderate, -conservative or -expert");
    had_err = -1;
  }
  else if (args->aggressive || args->state_of_truth)
  {
    args->hmin = 3UL;
    args->read_hmin = 1UL;
    args->altmax = (double) 1.00;
    args->refmin = (double) 0.00;
    args->mapqmin = 0UL;
    args->covmin = 1UL;
    args->clenmax = ULONG_MAX;
    args->allow_multiple = true;
  }
  else if (args->moderate)
  {
    args->hmin = 3UL;
    args->read_hmin = 1UL;
    args->altmax = (double) 0.99;
    args->refmin = (double) 0.00;
    args->mapqmin = 10UL;
    args->covmin = 1UL;
    args->clenmax = ULONG_MAX;
    args->allow_multiple = true;
  }
  else if (args->conservative)
  {
    args->hmin = 3UL;
    args->read_hmin = 2UL;
    args->altmax = (double) 0.80;
    args->refmin = (double) 0.10;
    args->mapqmin = 21UL;
    args->covmin = 1UL;
    args->clenmax = ULONG_MAX;
    args->allow_multiple = false;
  }

  if (gt_str_length(args->outprefix) == 0)
  {
    gt_error_set(err, "outprefix cannot be an empty string");
    had_err = -1;
  }
  return had_err;
}

static int gt_hop_runner(GT_UNUSED int argc, GT_UNUSED const char **argv,
    GT_UNUSED int parsed_args, void *tool_arguments, GtError *err)
{
  GtHopArguments *arguments = tool_arguments;
  GtEncseq *encseq;
  GtEncseqLoader *el;
  GtLogger *v_logger;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);
  v_logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  gt_logger_log(v_logger, "Correction parameters:");
  gt_logger_log(v_logger, "hmin = %lu", arguments->hmin);
  gt_logger_log(v_logger, "read-hmin = %lu", arguments->read_hmin);
  gt_logger_log(v_logger, "altmax = %.2f", arguments->altmax);
  gt_logger_log(v_logger, "refmin = %.2f", arguments->refmin);
  gt_logger_log(v_logger, "mapqmin = %lu", arguments->mapqmin);
  gt_logger_log(v_logger, "covmin = %lu", arguments->covmin);
  if (arguments->clenmax == ULONG_MAX)
    gt_logger_log(v_logger, "clenmax = unlimited");
  else
    gt_logger_log(v_logger, "clenmax = %lu", arguments->clenmax);
  gt_logger_log(v_logger, "allow-multiple = %s", arguments->allow_multiple
      ? "yes" : "no");
  if (gt_str_length(arguments->annotation) > 0)
    gt_logger_log(v_logger, "restrict to %s feature in annotation %s",
        gt_str_get(arguments->atype), gt_str_get(arguments->annotation));
  gt_assert(gt_str_length(arguments->encseqinput) > 0);
  el = gt_encseq_loader_new();
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
    GtSamfileEncseqMapping *sem = NULL;
    GtAlignedSegmentsPile *asp = NULL;
    GtFile *outfile = NULL;
    GtAlphabet *dna = gt_alphabet_new_dna();
    GtSeqIterator **readset_iters = NULL;
    GtStrArray **infiles = NULL;
    GtFile **outfiles = NULL;
    unsigned long nfiles = 0, i;
    hpp = gt_hpol_processor_new(encseq, arguments->hmin);
    if (gt_str_length(arguments->map) > 0)
    {
      if (!arguments->map_is_sam)
        sfi = gt_samfile_iterator_new_bam(gt_str_get(arguments->map), dna, err);
      else
        sfi = gt_samfile_iterator_new_sam(gt_str_get(arguments->map), dna,
            NULL, err);
      if (sfi == NULL)
        had_err = -1;
      sem = gt_samfile_encseq_mapping_new(sfi, encseq, err);
      if (sem == NULL)
        had_err = -1;
      if (!had_err)
      {
        asp = gt_aligned_segments_pile_new(sfi, sem);
        if (!arguments->rchk)
          gt_hpol_processor_enable_segments_hlen_adjustment(hpp, asp,
              arguments->read_hmin, arguments->altmax, arguments->refmin,
              arguments->mapqmin, arguments->covmin, arguments->allow_partial,
              arguments->allow_multiple, arguments->clenmax);
        else
          gt_hpol_processor_enable_aligned_segments_refregionscheck(hpp, asp);
        if (arguments->stats || arguments->state_of_truth)
          gt_hpol_processor_enable_statistics_output(hpp,
              arguments->state_of_truth, NULL);
        if (!had_err && gt_str_length(arguments->outfilename) > 0)
        {
          outfile = gt_file_new(gt_str_get(arguments->outfilename), "w", err);
          if (outfile == NULL)
            had_err = -1;
          else
            gt_hpol_processor_enable_direct_segments_output(hpp, outfile);
        }
        else if (!had_err)
        {
          nfiles = gt_str_array_size(arguments->readset);
          if (!had_err && nfiles > 0)
          {
            GtStr *outfn = gt_str_new();
            infiles = gt_calloc((size_t)nfiles, sizeof (*infiles));
            outfiles = gt_calloc((size_t)nfiles, sizeof (*outfiles));
            readset_iters = gt_malloc(sizeof (*readset_iters) * nfiles);
            for (i = 0; i < nfiles && !had_err; i++)
            {
              char *bn;
              gt_str_set(outfn, gt_str_get(arguments->outprefix));
              infiles[i] = gt_str_array_new();
              gt_str_array_add_cstr(infiles[i],
                  gt_str_array_get(arguments->readset, i));
              bn = gt_basename(gt_str_array_get(arguments->readset, i));
              gt_str_append_cstr(outfn, bn);
              outfiles[i] = gt_file_new(gt_str_get(outfn), "w", err);
              if (outfiles[i] == NULL)
                had_err = -1;
              readset_iters[i] =
                gt_seq_iterator_fastq_new(infiles[i], err);
              if (readset_iters[i] == NULL)
                had_err = -1;
              gt_seq_iterator_fastq_relax_check_of_quality_description(
                  (GtSeqIteratorFastQ*)readset_iters[i]);
              gt_free(bn);
            }
            if (!had_err)
              gt_hpol_processor_enable_sorted_segments_output(hpp, nfiles,
                  readset_iters, outfiles);
            gt_str_delete(outfn);
          }
        }
      }
    }
    if (!had_err && gt_str_length(arguments->annotation) > 0)
    {
      spc = gt_seqpos_classifier_new(gt_str_get(arguments->annotation),
          gt_str_get(arguments->atype));
      gt_hpol_processor_restrict_to_feature_type(hpp, spc);
    }
    if (!had_err)
      had_err = gt_hpol_processor_run(hpp, v_logger, err);
    gt_aligned_segments_pile_delete(asp);
    gt_samfile_iterator_delete(sfi);
    gt_samfile_encseq_mapping_delete(sem);
    gt_seqpos_classifier_delete(spc);
    gt_file_delete(outfile);
    gt_hpol_processor_delete(hpp);
    for (i = 0; i < nfiles; i++)
    {
      gt_assert(infiles != NULL);
      gt_assert(readset_iters != NULL);
      gt_assert(outfiles != NULL);
      gt_str_array_delete(infiles[i]);
      gt_seq_iterator_delete(readset_iters[i]);
      gt_file_delete(outfiles[i]);
    }
    gt_free(infiles);
    gt_free(outfiles);
    gt_free(readset_iters);
    gt_alphabet_delete(dna);
  }
  gt_logger_delete(v_logger);
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  return had_err;
}

GtTool* gt_hop(void)
{
  return gt_tool_new(gt_hop_arguments_new,
                  gt_hop_arguments_delete,
                  gt_hop_option_parser_new,
                  gt_hop_arguments_check,
                  gt_hop_runner);
}
