/*
  Copyright (c) 2008-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2013 Center for Bioinformatics, University of Hamburg

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

#include "core/bioseq.h"
#include "core/encseq.h"
#include "core/fileutils_api.h"
#include "core/log.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/range.h"
#include "core/safearith.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/region_mapping.h"
#include "extended/seqid2file.h"
#include "extended/visitor_stream.h"
#include "ltr/gt_ltrdigest.h"
#include "ltr/ltr_input_check_visitor.h"
#include "ltr/ltrdigest_def.h"
#include "ltr/ltrdigest_file_out_stream.h"
#include "ltr/ltrdigest_pbs_visitor.h"
#include "ltr/ltrdigest_pdom_visitor.h"
#include "ltr/ltrdigest_ppt_visitor.h"
#include "ltr/ltrdigest_strand_assign_visitor.h"
#include "ltr/pdom_model_set.h"

typedef struct GtLTRdigestOptions {
  GtStr *trna_lib, *prefix, *cutoffs;
  bool verbose,
       write_alignments,
       write_aaseqs,
       output_all_chains,
       print_metadata,
       force_recreate;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
  GtStrArray *hmm_files;
  GtSeqid2FileInfo *s2fi;
  GtPdomCutoff cutoff;
  double evalue_cutoff;
  GtUword nthreads;
  unsigned int chain_max_gap_length,
               seqnamelen;
  GtRange ppt_len, ubox_len;
  double ppt_pyrimidine_prob,
         ppt_purine_prob,
         bkg_a_prob,
         bkg_g_prob,
         bkg_t_prob,
         bkg_c_prob,
         ubox_u_prob;
  unsigned int ppt_radius,
               max_ubox_dist;
  unsigned int pbs_radius,
               max_edist;
  GtRange alilen,
          offsetlen,
          trnaoffsetlen;
  int ali_score_match,
      ali_score_mismatch,
      ali_score_insertion,
      ali_score_deletion;
  GtBioseq *trna_lib_bs;
} GtLTRdigestOptions;

static void* gt_ltrdigest_arguments_new(void)
{
  GtLTRdigestOptions *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  memset(arguments, 0, sizeof *arguments);
  arguments->trna_lib = gt_str_new();
  arguments->prefix = gt_str_new();
  arguments->cutoffs = gt_str_new();
  arguments->ofi = gt_output_file_info_new();
  arguments->hmm_files = gt_str_array_new();
  arguments->s2fi = gt_seqid2file_info_new();
  return arguments;
}

static void gt_ltrdigest_arguments_delete(void *tool_arguments)
{
  GtLTRdigestOptions *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->trna_lib);
  gt_str_delete(arguments->prefix);
  gt_str_delete(arguments->cutoffs);
  gt_str_array_delete(arguments->hmm_files);
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_free(arguments);
}

static GtOptionParser* gt_ltrdigest_option_parser_new(void *tool_arguments)
{
  GtLTRdigestOptions *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *o, *ot, *oto;
  GtOption *oh, *oc, *oeval;
  static const char *cutoffs[] = {"NONE", "GA", "TC", NULL};
  static GtRange pptlen_defaults           = { 8UL, 30UL},
                 uboxlen_defaults          = { 3UL, 30UL},
                 pbsalilen_defaults        = {11UL, 30UL},
                 pbsoffsetlen_defaults     = { 0UL,  5UL},
                 pbstrnaoffsetlen_defaults = { 0UL,  5UL};
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] gff3_file",
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

  o = gt_option_new_bool("metadata",
                         "output metadata (run conditions) to separate file",
                         &arguments->print_metadata,
                         true);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, oto);

  o = gt_option_new_uint("seqnamelen",
                         "set maximal length of sequence names in FASTA headers"
                         " (e.g. for clustalw or similar tools)",
                         &arguments->seqnamelen,
                         20U);
  gt_option_parser_add_option(op, o);

  /* PPT search options */

  o = gt_option_new_range("pptlen",
                          "required PPT length range",
                          &arguments->ppt_len,
                          &pptlen_defaults);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_range("uboxlen",
                          "required U-box length range",
                          &arguments->ubox_len,
                          &uboxlen_defaults);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_uint("uboxdist",
                         "allowed U-box distance range from PPT",
                         &arguments->max_ubox_dist, 0);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_uint("pptradius",
                         "radius around beginning of 3' LTR "
                         "to search for PPT",
                         &arguments->ppt_radius,
                         30U);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_probability("pptrprob",
                                "purine emission probability inside PPT",
                                &arguments->ppt_purine_prob,
                                PPT_PURINE_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("pptyprob",
                                "pyrimidine emission probability inside PPT",
                                &arguments->ppt_pyrimidine_prob,
                                PPT_PYRIMIDINE_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("pptgprob",
                                "background G emission probability outside PPT",
                                &arguments->bkg_g_prob,
                                BKG_G_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("pptcprob",
                                "background C emission probability outside PPT",
                                &arguments->bkg_c_prob,
                                BKG_C_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("pptaprob",
                                "background A emission probability outside PPT",
                                &arguments->bkg_a_prob,
                                BKG_A_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("ppttprob",
                                "background T emission probability outside PPT",
                                &arguments->bkg_t_prob,
                                BKG_T_PROB);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  o = gt_option_new_probability("pptuprob",
                                "U/T emission probability inside U-box",
                                &arguments->ubox_u_prob,
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
                          &arguments->alilen,
                          &pbsalilen_defaults);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

  o = gt_option_new_range("pbsoffset",
                          "allowed PBS offset from LTR boundary range",
                          &arguments->offsetlen,
                          &pbsoffsetlen_defaults);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

  o = gt_option_new_range("pbstrnaoffset",
                          "allowed PBS/tRNA 3' end alignment offset range",
                          &arguments->trnaoffsetlen,
                          &pbstrnaoffsetlen_defaults);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

  o = gt_option_new_uint("pbsmaxedist",
                         "maximal allowed PBS/tRNA alignment unit "
                         "edit distance",
                         &arguments->max_edist,
                         1U);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

  o = gt_option_new_uint("pbsradius",
                         "radius around end of 5' LTR "
                         "to search for PBS",
                         &arguments->pbs_radius,
                         30U);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

 /* Protein domain search options */

  oh = gt_option_new_filename_array("hmms",
                                    "profile HMM models for domain detection "
                                    "(separate by spaces, finish with --) in "
                                    "HMMER3 format\n"
                                    "Omit this option to disable pHMM search.",
                                    arguments->hmm_files);
  gt_option_parser_add_option(op, oh);

  oeval = gt_option_new_probability("pdomevalcutoff",
                                    "global E-value cutoff for pHMM search\n"
                                    "default 1E-6",
                                    &arguments->evalue_cutoff,
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
                         &arguments->write_alignments,
                         false);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, oh);
  gt_option_imply(o, oto);

  o = gt_option_new_bool("aaout",
                         "output amino acid sequences for protein domain "
                         "hits",
                         &arguments->write_aaseqs,
                         false);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, oh);
  gt_option_imply(o, oto);

  o = gt_option_new_bool("allchains",
                           "output features from all chains and unchained "
                           "features, labeled with chain numbers",
                           &arguments->output_all_chains,
                           false);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, oh);

  o = gt_option_new_uint("maxgaplen",
                         "maximal allowed gap size between fragments (in amino "
                         "acids) when chaining pHMM hits for a protein domain",
                         &arguments->chain_max_gap_length,
                         50U);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
  gt_option_imply(o, oh);

  o = gt_option_new_uword("threads",
                          "DEPRECATED, only included for compatibility reasons!"
                          " Use the -j parameter of the 'gt' call instead.",
                          &arguments->nthreads,
                          0);
  gt_option_parser_add_option(op, o);
  gt_option_is_development_option(o);

  o = gt_option_new_bool("force_recreate",
                         "force recreation of hmmpressed profiles",
                         &arguments->force_recreate,
                         false);
  gt_option_parser_add_option(op, o);

  /* Extended PBS options */

  o = gt_option_new_int("pbsmatchscore",
                        "match score for PBS/tRNA alignments",
                        &arguments->ali_score_match,
                        5);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
  gt_option_imply(o, ot);

  o = gt_option_new_int("pbsmismatchscore",
                        "mismatch score for PBS/tRNA alignments",
                        &arguments->ali_score_mismatch,
                        -10);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
  gt_option_imply(o, ot);

  o = gt_option_new_int("pbsinsertionscore",
                        "insertion score for PBS/tRNA alignments",
                        &arguments->ali_score_insertion,
                        -20);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
  gt_option_imply(o, ot);

  o = gt_option_new_int("pbsdeletionscore",
                        "deletion score for PBS/tRNA alignments",
                        &arguments->ali_score_deletion,
                        -20);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
  gt_option_imply(o, ot);

  /* verbosity */

  o = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, o);

  /* output file options */

  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  /* region mapping and sequence source options */

  gt_seqid2file_register_options_ext(op, arguments->s2fi, false, false);

  return op;
}

int gt_ltrdigest_arguments_check(GT_UNUSED int rest_argc, void *tool_arguments,
                                 GtError* err)
{
  GtLTRdigestOptions *arguments = tool_arguments;
  int had_err  = 0;

  if (arguments->nthreads > 0) {
    gt_warning("The '-threads' option is deprecated. Please use the '-j'"
               "option of the 'gt' call instead, e.g.:\n"
               "  gt -j "GT_WU" ltrdigest ...", arguments->nthreads);
  }

  /* -trnas */
  if (!had_err && arguments->trna_lib && gt_str_length(arguments->trna_lib) > 0)
  {
    if (!gt_file_exists(gt_str_get(arguments->trna_lib)))
    {
      gt_error_set(err, "File '%s' does not exist!",
                        gt_str_get(arguments->trna_lib));
      had_err = -1;
    }
  }
  return had_err;
}

static int gt_ltrdigest_runner(GT_UNUSED int argc, const char **argv,
                               int parsed_args, void *tool_arguments,
                               GtError *err)
{
  GtLTRdigestOptions *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream  = NULL,
               *check_stream    = NULL,
               *gff3_out_stream = NULL,
               *pdom_stream     = NULL,
               *ppt_stream      = NULL,
               *pbs_stream      = NULL,
               *tab_out_stream  = NULL,
               *sa_stream       = NULL,
               *last_stream     = NULL;
  int had_err      = 0,
      tests_to_run = 0,
      arg = parsed_args;
  GtRegionMapping *rmap = NULL;
  GtPdomModelSet *ms = NULL;
  gt_error_check(err);
  gt_assert(arguments);

  /* determine and open sequence source */
  if (gt_seqid2file_option_used(arguments->s2fi)) {
    /* create region mapping */
    rmap = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
    if (!rmap)
      had_err = -1;
  } else {
    GtEncseqLoader *el;
    GtEncseq *encseq;
    /* no new-style sequence source option given, fall back to legacy syntax */
    if (argc < 3) {
      gt_error_set(err, "missing mandatory argument(s)");
      had_err = -1;
    }
    if (!had_err) {
      el = gt_encseq_loader_new();
      gt_encseq_loader_disable_autosupport(el);
      gt_encseq_loader_require_md5_support(el);
      gt_encseq_loader_require_description_support(el);
      encseq = gt_encseq_loader_load(el, argv[argc-1], err);
      /* XXX: clip off terminal argument */
      gt_free((char*) argv[argc-1]);
      argv[argc-1] = NULL;
      argc--;
      gt_encseq_loader_delete(el);
      if (!encseq)
        had_err = -1;
      else {
        rmap = gt_region_mapping_new_encseq_seqno(encseq);
        gt_encseq_delete(encseq);
      }
    }
  }
  gt_assert(had_err || rmap);

  /* Always search for PPT. */
  tests_to_run |= GT_LTRDIGEST_RUN_PPT;

  /* Open tRNA library if given. */
  if (!had_err && arguments->trna_lib
        && gt_str_length(arguments->trna_lib) > 0)
  {
    tests_to_run |= GT_LTRDIGEST_RUN_PBS;
    arguments->trna_lib_bs = gt_bioseq_new(gt_str_get(arguments->trna_lib),
                                           err);
    if (gt_error_is_set(err))
      had_err = -1;
  }

  /* Set HMMER cutoffs. */
  if (!had_err && gt_str_array_size(arguments->hmm_files) > 0)
  {
    tests_to_run |= GT_LTRDIGEST_RUN_PDOM;
    if (!strcmp(gt_str_get(arguments->cutoffs), "GA")) {
      arguments->cutoff = GT_PHMM_CUTOFF_GA;
    } else if (!strcmp(gt_str_get(arguments->cutoffs), "TC")) {
      arguments->cutoff = GT_PHMM_CUTOFF_TC;
    } else if (!strcmp(gt_str_get(arguments->cutoffs), "NONE")) {
      arguments->cutoff = GT_PHMM_CUTOFF_NONE;
    } else {
      gt_error_set(err, "invalid cutoff setting!");
      had_err = -1;
    }
  }

  if (!had_err) {
    last_stream = gff3_in_stream  = gt_gff3_in_stream_new_sorted(argv[arg]);
  }

  if (!had_err) {
    GtNodeVisitor *check_v;
    check_v = gt_ltr_input_check_visitor_new();
    last_stream = check_stream = gt_visitor_stream_new(last_stream, check_v);
  }

  if (!had_err && gt_str_array_size(arguments->hmm_files) > 0) {
    GtNodeVisitor *pdom_v;
    ms = gt_pdom_model_set_new(arguments->hmm_files, arguments->force_recreate,
                               err);
    if (ms != NULL) {
      pdom_v = gt_ltrdigest_pdom_visitor_new(ms, arguments->evalue_cutoff,
                                             arguments->chain_max_gap_length,
                                             arguments->cutoff, rmap, err);
      if (pdom_v == NULL)
        had_err = -1;
      if (!had_err) {
        gt_ltrdigest_pdom_visitor_set_source_tag((GtLTRdigestPdomVisitor*)
                                                                        pdom_v,
                                                 GT_LTRDIGEST_TAG);
        if (arguments->output_all_chains)
          gt_ltrdigest_pdom_visitor_output_all_chains((GtLTRdigestPdomVisitor*)
                                                                        pdom_v);
        last_stream = pdom_stream = gt_visitor_stream_new(last_stream, pdom_v);
      }
    } else had_err = -1;
  }

  if (!had_err && arguments->trna_lib_bs) {
    GtNodeVisitor *pbs_v;
    pbs_v = gt_ltrdigest_pbs_visitor_new(rmap, arguments->pbs_radius,
                                         arguments->max_edist,
                                         arguments->alilen,
                                         arguments->offsetlen,
                                         arguments->trnaoffsetlen,
                                         arguments->ali_score_match,
                                         arguments->ali_score_mismatch,
                                         arguments->ali_score_insertion,
                                         arguments->ali_score_deletion,
                                         arguments->trna_lib_bs, err);
    if (pbs_v != NULL)
      last_stream = pbs_stream = gt_visitor_stream_new(last_stream, pbs_v);
    else
      had_err = -1;
  }

  if (!had_err) {
    GtNodeVisitor *ppt_v;
    ppt_v = gt_ltrdigest_ppt_visitor_new(rmap, arguments->ppt_len,
                                         arguments->ubox_len,
                                         arguments->ppt_pyrimidine_prob,
                                         arguments->ppt_purine_prob,
                                         arguments->bkg_a_prob,
                                         arguments->bkg_g_prob,
                                         arguments->bkg_t_prob,
                                         arguments->bkg_c_prob,
                                         arguments->ubox_u_prob,
                                         arguments->ppt_radius,
                                         arguments->max_ubox_dist, err);
    if (ppt_v != NULL)
      last_stream = ppt_stream = gt_visitor_stream_new(last_stream, ppt_v);
    else
      had_err = -1;
  }

  if (!had_err) {
    GtNodeVisitor *sa_v;
    sa_v = gt_ltrdigest_strand_assign_visitor_new();
    gt_assert(sa_v);
    last_stream = sa_stream = gt_visitor_stream_new(last_stream, sa_v);
  }

  if (!had_err)
  {
    /* attach tabular output stream, if requested */
    if (gt_str_length(arguments->prefix) > 0)
    {
      last_stream = tab_out_stream = gt_ltrdigest_file_out_stream_new(
                                                  last_stream,
                                                  tests_to_run,
                                                  rmap,
                                                  gt_str_get(arguments->prefix),
                                                  arguments->seqnamelen,
                                                  err);
      if (!tab_out_stream)
        had_err = -1;
      if (!had_err && arguments->print_metadata)
      {
        had_err = gt_ltrdigest_file_out_stream_write_metadata(
                                           (GtLTRdigestFileOutStream*)
                                                                 tab_out_stream,
                                           tests_to_run,
                                           gt_str_get(arguments->trna_lib),
                                           argv[arg],
                                           arguments->ppt_len,
                                           arguments->ubox_len,
                                           arguments->ppt_radius,
                                           arguments->alilen,
                                           arguments->max_edist,
                                           arguments->offsetlen,
                                           arguments->trnaoffsetlen,
                                           arguments->pbs_radius,
                                           arguments->hmm_files,
                                           arguments->chain_max_gap_length,
                                           arguments->evalue_cutoff,
                                           err);
      }
      if (!had_err)
      {
        if (arguments->write_alignments)
          gt_ltrdigest_file_out_stream_enable_pdom_alignment_output(
                                                                tab_out_stream);
        if (arguments->write_aaseqs)
          gt_ltrdigest_file_out_stream_enable_aa_sequence_output(
                                                                tab_out_stream);
      }
    }

    if (!had_err) {
      last_stream = gff3_out_stream = gt_gff3_out_stream_new(last_stream,
                                                             arguments->outfp);
      gt_assert(last_stream);

      /* pull the features through the stream and free them afterwards */
      had_err = gt_node_stream_pull(last_stream, err);
    }
  }

  gt_pdom_model_set_delete(ms);
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(ppt_stream);
  gt_node_stream_delete(pbs_stream);
  gt_node_stream_delete(sa_stream);
  gt_node_stream_delete(pdom_stream);
  gt_node_stream_delete(tab_out_stream);
  gt_node_stream_delete(check_stream);
  gt_node_stream_delete(gff3_in_stream);
  gt_bioseq_delete(arguments->trna_lib_bs);
  gt_region_mapping_delete(rmap);

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
