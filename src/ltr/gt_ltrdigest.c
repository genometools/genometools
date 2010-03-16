/*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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
#include "core/fileutils_api.h"
#include "core/log.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/option.h"
#include "core/outputfile.h"
#include "core/safearith.h"
#include "core/unused_api.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "ltr/gt_ltrdigest.h"
#include "ltr/ltrdigest_def.h"
#include "ltr/ltrdigest_stream.h"
#include "ltr/ltrfileout_stream.h"
#include "match/encodedsequence.h"

typedef struct GtLTRdigestOptions {
  GtPBSOptions  pbs_opts;
  GtPPTOptions  ppt_opts;
#ifdef HAVE_HMMER
  GtPdomOptions pdom_opts;
#endif
  GtStr *trna_lib, *prefix, *cutoffs;
  bool verbose;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
  unsigned int seqnamelen;
} GtLTRdigestOptions;

static void* gt_ltrdigest_arguments_new(void)
{
  GtLTRdigestOptions *arguments = gt_calloc(1, sizeof *arguments);
  memset(arguments, 0, sizeof *arguments);
#ifdef HAVE_HMMER
  arguments->pdom_opts.hmm_files = gt_str_array_new();
#endif
  arguments->trna_lib = gt_str_new();
  arguments->prefix = gt_str_new();
  arguments->cutoffs = gt_str_new();
  arguments->ofi = gt_outputfileinfo_new();
  return arguments;
}

static void gt_ltrdigest_arguments_delete(void *tool_arguments)
{
  GtLTRdigestOptions *arguments = tool_arguments;
  if (!arguments) return;
#ifdef HAVE_HMMER
  gt_str_array_delete(arguments->pdom_opts.hmm_files);
#endif
  gt_str_delete(arguments->trna_lib);
  gt_str_delete(arguments->prefix);
  gt_str_delete(arguments->cutoffs);
  gt_file_delete(arguments->outfp);
  gt_outputfileinfo_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_ltrdigest_option_parser_new(void *tool_arguments)
{
  GtLTRdigestOptions *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *o, *ot, *oto;
#ifdef HAVE_HMMER
  GtOption *oh, *oc, *oeval;
  static const char *cutoffs[] = {"NONE", "GA", "TC", NULL};
#endif
  static GtRange pptlen_defaults           = { 8, 30},
                 uboxlen_defaults          = { 3, 30},
                 pbsalilen_defaults        = {11, 30},
                 pbsoffsetlen_defaults     = { 0,  5},
                 pbstrnaoffsetlen_defaults = { 0,  5};
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] gff3_file indexname",
                            "Identifies and annotates sequence features in LTR "
                            "retrotransposon candidates.");

  /* Output files */

  oto = gt_option_new_string("outfileprefix",
                             "prefix for output files (e.g. 'foo' will create "
                             "files called 'foo_*.csv' and 'foo_*.fas')\n"
                             "Omit this option for GFF3 output only.",
                             arguments->prefix,
                             NULL);
  gt_option_parser_add_option(op, oto);
  gt_option_hide_default(oto);

  o = gt_option_new_uint("seqnamelen",
                         "set maximal length of sequence names in FASTA headers"
                         " (e.g. for clustalw or similar tools)",
                         &arguments->seqnamelen,
                         20);
  gt_option_parser_add_option(op, o);

  /* PPT search options */

  o = gt_option_new_range("pptlen",
                          "required PPT length range",
                          &arguments->ppt_opts.ppt_len,
                          &pptlen_defaults);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_range("uboxlen",
                          "required U-box length range",
                          &arguments->ppt_opts.ubox_len,
                          &uboxlen_defaults);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_uint("pptradius",
                         "radius around beginning of 3' LTR "
                         "to search for PPT",
                         &arguments->ppt_opts.radius,
                         30);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_probability("pptrprob",
                                "purine emission probability inside PPT",
                                &arguments->ppt_opts.ppt_purine_prob,
                                PPT_PURINE_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("pptyprob",
                                "pyrimidine emission probability inside PPT",
                                &arguments->ppt_opts.ppt_pyrimidine_prob,
                                PPT_PYRIMIDINE_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("pptgprob",
                                "background G emission probability outside PPT",
                                &arguments->ppt_opts.bkg_g_prob,
                                BKG_G_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("pptcprob",
                                "background C emission probability outside PPT",
                                &arguments->ppt_opts.bkg_c_prob,
                                BKG_C_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("pptaprob",
                                "background A emission probability outside PPT",
                                &arguments->ppt_opts.bkg_a_prob,
                                BKG_A_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("ppttprob",
                                "background T emission probability outside PPT",
                                &arguments->ppt_opts.bkg_t_prob,
                                BKG_T_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("pptuprob",
                                "U/T emission probability inside U-box",
                                &arguments->ppt_opts.ubox_u_prob,
                                UBOX_U_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  /* PBS search options */

  ot = gt_option_new_filename("trnas",
                              "tRNA library in multiple FASTA format for PBS "
                              "detection\n"
                              "Omit this option to disable PBS search.",
                              arguments->trna_lib);
  gt_option_parser_add_option(op, ot);
  gt_option_hide_default(ot);

  o = gt_option_new_range("pbsalilen",
                          "required PBS/tRNA alignment length range",
                          &arguments->pbs_opts.alilen,
                          &pbsalilen_defaults);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

  o = gt_option_new_range("pbsoffset",
                          "allowed PBS offset from LTR boundary range",
                          &arguments->pbs_opts.offsetlen,
                          &pbsoffsetlen_defaults);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

  o = gt_option_new_range("pbstrnaoffset",
                          "allowed PBS/tRNA 3' end alignment offset range",
                          &arguments->pbs_opts.trnaoffsetlen,
                          &pbstrnaoffsetlen_defaults);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

  o = gt_option_new_uint("pbsmaxedist",
                         "maximal allowed PBS/tRNA alignment unit "
                         "edit distance",
                         &arguments->pbs_opts.max_edist,
                         1);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

  o = gt_option_new_uint("pbsradius",
                         "radius around end of 5' LTR "
                         "to search for PBS",
                         &arguments->pbs_opts.radius,
                         30);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

 /* Protein domain search options */

#ifdef HAVE_HMMER
  oh = gt_option_new_filenamearray("hmms",
                                   "profile HMM models for domain detection "
                                   "(separate by spaces, finish with --) in "
                                   "HMMER2 format\n"
                                   "Omit this option to disable pHMM search.",
                                   arguments->pdom_opts.hmm_files);
  gt_option_parser_add_option(op, oh);

  oeval = gt_option_new_probability("pdomevalcutoff",
                                    "global E-value cutoff for pHMM search\n"
                                    "default 1E-6",
                                    &arguments->pdom_opts.evalue_cutoff,
                                    0.000001);
  gt_option_parser_add_option(op, oeval);
  gt_option_is_extended_option(oeval);
  gt_option_hide_default(oeval);
  gt_option_imply(oeval, oh);

  oc = gt_option_new_choice("pdomcutoff", "model-specific score cutoff\n"
                                       "choose from TC (trusted cutoff) | "
                                       "GA (gathering cutoff) | "
                                       "NONE (no cutoffs)",
                             arguments->cutoffs, cutoffs[0], cutoffs);
  gt_option_parser_add_option(op, oc);
  gt_option_is_extended_option(oeval);
  gt_option_imply(oeval, oh);

  o = gt_option_new_bool("aliout",
                           "output pHMM to amino acid sequence alignments",
                           &arguments->pdom_opts.write_alignments,
                           false);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, oh);
  gt_option_imply(o, oto);

  o = gt_option_new_bool("aaout",
                           "output amino acid sequences for protein domain "
                           "hits",
                           &arguments->pdom_opts.write_aaseqs,
                           false);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, oh);
  gt_option_imply(o, oto);

  o = gt_option_new_uint_min("threads",
                             "number of concurrent worker threads to use in "
                             "pHMM scanning",
                             &arguments->pdom_opts.nof_threads,
                             2, 1);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, oh);

  o = gt_option_new_uint("maxgaplen",
                         "maximal allowed gap size between fragments (in amino "
                         "acids) when chaining pHMM hits for a protein domain",
                         &arguments->pdom_opts.chain_max_gap_length,
                         50);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
  gt_option_imply(o, oh);
#endif

  /* Extended PBS options */

  o = gt_option_new_int("pbsmatchscore",
                        "match score for PBS/tRNA alignments",
                        &arguments->pbs_opts.ali_score_match,
                        5);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
  gt_option_imply(o, ot);

  o = gt_option_new_int("pbsmismatchscore",
                        "mismatch score for PBS/tRNA alignments",
                        &arguments->pbs_opts.ali_score_mismatch,
                        -10);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
  gt_option_imply(o, ot);

  o = gt_option_new_int("pbsinsertionscore",
                        "insertion score for PBS/tRNA alignments",
                        &arguments->pbs_opts.ali_score_insertion,
                        -20);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
  gt_option_imply(o, ot);

  o = gt_option_new_int("pbsdeletionscore",
                        "deletion score for PBS/tRNA alignments",
                        &arguments->pbs_opts.ali_score_deletion,
                        -20);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
  gt_option_imply(o, ot);

  /* verbosity */

  o = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, o);

  /* output file options */

  gt_outputfile_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_min_max_args(op, 2, 2);

  return op;
}

int gt_ltrdigest_arguments_check(GT_UNUSED int rest_argc, void *tool_arguments,
                                 GtError* err)
{
  GtLTRdigestOptions *arguments = tool_arguments;
  int had_err  = 0;

  /* -trnas */
  if (arguments->trna_lib
        && gt_str_length(arguments->trna_lib) > 0)
  {
    if (!gt_file_exists(gt_str_get(arguments->trna_lib)))
    {
      gt_error_set(err, "File '%s' does not exist!",
                        gt_str_get(arguments->trna_lib));
      had_err = -1;
    }
  }

  if (!had_err)
  {
    GtHMM *hmm;
    GtAlphabet *alpha;
    alpha = gt_alphabet_new_dna();
    hmm = gt_ppt_hmm_new(alpha, &arguments->ppt_opts);
    if (!hmm)
    {
      gt_error_set(err, "PPT HMM parameters are not valid!");
      had_err = -1;
    }
    else
      gt_hmm_delete(hmm);
    gt_alphabet_delete(alpha);
  }

  return had_err;
}

static int gt_ltrdigest_runner(GT_UNUSED int argc, const char **argv,
                               int parsed_args, void *tool_arguments,
                               GtError *err)
{
  GtLTRdigestOptions *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream   = NULL,
               *gff3_out_stream  = NULL,
               *ltrdigest_stream = NULL,
               *tab_out_stream   = NULL,
               *last_stream      = NULL;
  int had_err      = 0,
      tests_to_run = 0,
      arg = parsed_args;
  GtStr *indexname = gt_str_new_cstr(argv[arg+1]);
  GtLogger *logger = gt_logger_new(arguments->verbose,
                                   GT_LOGGER_DEFLT_PREFIX, stdout);
  gt_error_check(err);
  gt_assert(arguments);

  /* Open sequence file */
  GtEncodedsequence *encseq = gt_encodedsequence_new_from_index(true,
                                               indexname,
                                               true,
                                               true,
                                               true,
                                               true,
                                               logger,
                                               err);
  if (gt_error_is_set(err))
    had_err = -1;

  /* Always search for PPT. */
  tests_to_run |= GT_LTRDIGEST_RUN_PPT;

  /* Open tRNA library if given. */
  if (!had_err && arguments->trna_lib
        && gt_str_length(arguments->trna_lib) > 0)
  {
    tests_to_run |= GT_LTRDIGEST_RUN_PBS;
   arguments->pbs_opts.trna_lib = gt_bioseq_new(gt_str_get(arguments->trna_lib),
                                                 err);
    if (gt_error_is_set(err))
      had_err = -1;
  }

#ifdef HAVE_HMMER
  /* Open HMMER files if given. */
  if (!had_err && gt_str_array_size(arguments->pdom_opts.hmm_files) > 0)
  {
    tests_to_run |= GT_LTRDIGEST_RUN_PDOM;
    if (!strcmp(gt_str_get(arguments->cutoffs), "GA")) {
      arguments->pdom_opts.cutoff = GT_PHMM_CUTOFF_GA;
    } else if (!strcmp(gt_str_get(arguments->cutoffs), "TC")) {
      arguments->pdom_opts.cutoff = GT_PHMM_CUTOFF_TC;
    } else if (!strcmp(gt_str_get(arguments->cutoffs), "NONE")) {
      arguments->pdom_opts.cutoff = GT_PHMM_CUTOFF_NONE;
    } else {
      gt_error_set(err, "invalid cutoff setting!");
      had_err = -1;
    }
  }
#endif

  if (!had_err)
  {
    /* set up stream flow
     * ------------------*/
    last_stream = gff3_in_stream  = gt_gff3_in_stream_new_sorted(argv[arg]);

    last_stream = ltrdigest_stream = gt_ltrdigest_stream_new(last_stream,
                                                  tests_to_run,
                                                  encseq,
                                                  &arguments->pbs_opts,
                                                  &arguments->ppt_opts,
#ifdef HAVE_HMMER
                                                  &arguments->pdom_opts,
#endif
                                                  err);
    if (!ltrdigest_stream)
      had_err = -1;
  }

  if (!had_err)
  {
    /* attach tabular output stream, if requested */
    if (gt_str_length(arguments->prefix) > 0)
    {
      last_stream = tab_out_stream = gt_ltr_fileout_stream_new(last_stream,
                                              tests_to_run,
                                              encseq,
                                              gt_str_get(arguments->prefix),
                                              &arguments->ppt_opts,
                                              &arguments->pbs_opts,
#ifdef HAVE_HMMER
                                              &arguments->pdom_opts,
#endif
                                              gt_str_get(arguments->trna_lib),
                                              argv[arg+1],
                                              argv[arg],
                                              arguments->seqnamelen,
                                              err);
#ifdef HAVE_HMMER
    if (&arguments->pdom_opts.write_alignments)
      gt_ltr_fileout_stream_enable_pdom_alignment_output(tab_out_stream);
    if (&arguments->pdom_opts.write_aaseqs)
      gt_ltr_fileout_stream_enable_aa_sequence_output(tab_out_stream);
#endif
    }

    last_stream = gff3_out_stream = gt_gff3_out_stream_new(last_stream,
                                                           arguments->outfp);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(last_stream, err);
  }

  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(ltrdigest_stream);
  if (tab_out_stream)
    gt_node_stream_delete(tab_out_stream);
  gt_node_stream_delete(gff3_in_stream);

  gt_str_delete(indexname);
  gt_encodedsequence_delete(encseq);
  encseq = NULL;
  gt_bioseq_delete(arguments->pbs_opts.trna_lib);
  gt_logger_delete(logger);

  return had_err;
}

GtTool* gt_ltrdigest(void)
{
  return gt_tool_new(gt_ltrdigest_arguments_new,
                     gt_ltrdigest_arguments_delete,
                     gt_ltrdigest_option_parser_new,
                     gt_ltrdigest_arguments_check,
                     gt_ltrdigest_runner);
}
