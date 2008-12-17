/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include <ctype.h>
#include <string.h>
#include "core/bioseq.h"
#include "core/fileutils.h"
#include "core/ma.h"
#include "core/option.h"
#include "core/outputfile.h"
#include "core/safearith.h"
#include "core/unused_api.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream.h"
#include "extended/gtdatahelp.h"
#include "extended/seqid2file.h"
#include "ltr/gt_ltrdigest.h"
#include "ltr/ltrdigest_def.h"
#include "ltr/ltrdigest_stream.h"
#include "ltr/ltrfileout_stream.h"

typedef struct GtLTRdigestOptions {
  GtPBSOptions  pbs_opts;
  GtPPTOptions  ppt_opts;
#ifdef HAVE_HMMER
  GtPdomOptions pdom_opts;
#endif
  GtStr *trna_lib,
        *prefix,
        *regionmapping,
        *seqfile;
  bool verbose;
  GtOutputFileInfo *ofi;
  GtGenFile *outfp;
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
  arguments->seqfile = gt_str_new();
  arguments->regionmapping = gt_str_new();
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
  gt_str_delete(arguments->regionmapping);
  gt_str_delete(arguments->seqfile);
  gt_str_delete(arguments->trna_lib);
  gt_str_delete(arguments->prefix);
  gt_genfile_close(arguments->outfp);
  gt_outputfileinfo_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_ltrdigest_option_parser_new(void *tool_arguments)
{
  GtLTRdigestOptions *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *o, *ot, *oto;
#ifdef HAVE_HMMER
  GtOption *oh;
#endif
  static GtRange pptlen_defaults           = { 8, 30},
                 uboxlen_defaults          = { 3, 30},
                 pbsalilen_defaults        = {11, 30},
                 pbsoffsetlen_defaults     = { 0,  5},
                 pbstrnaoffsetlen_defaults = { 0,  5};
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [gff3_file]",
                         "Discovers and annotates sequence features in LTR "
                         "retrotransposon candidates.");

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
                      "maximal allowed PBS/tRNA alignment unit edit distance",
                      &arguments->pbs_opts.max_edist,
                      1);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

  o = gt_option_new_uint("pbsradius",
                      "radius around end of 5' LTR "
                      " to search for PBS",
                      &arguments->pbs_opts.radius,
                      30);
  gt_option_parser_add_option(op, o);
  gt_option_imply(o, ot);

 /* Protein domain search options */
#ifdef HAVE_HMMER
  oh = gt_option_new_filenamearray("hmms",
                               "profile HMM models for domain detection "
                               "(separate by spaces, finish with --) in HMMER"
                               "2 format\n"
                               "Omit this option to disable pHMM search.",
                               arguments->pdom_opts.hmm_files);
  gt_option_parser_add_option(op, oh);

  o = gt_option_new_probability("pdomevalcutoff",
                             "E-value cutoff for pHMM search",
                             &arguments->pdom_opts.evalue_cutoff,
                             0.000001);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
  gt_option_imply(o, oh);

  o = gt_option_new_uint_min("threads",
                          "number of concurrent worker threads to use in "
                          "pHMM scanning",
                          &arguments->pdom_opts.nof_threads,
                          2, 1);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_uint("maxgaplen",
                      "maximal allowed gap size between fragments (in amino "
                      "acids) when chaining pHMM hits for a protein domain",
                      &arguments->pdom_opts.chain_max_gap_length,
                      50);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);
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
                      "set maximal length of sequence names in FASTA headers "
                      "(e.g. for clustalw or similar tools)",
                      &arguments->seqnamelen,
                      20);
  gt_option_parser_add_option(op, o);

  /* Add -seqfile and -regionmapping arguments for consistency */
  gt_seqid2file_options(op, arguments->seqfile, arguments->regionmapping);

  /* verbosity */
  o = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, o);

  /* output file options */
  gt_outputfile_register_options(op, &arguments->outfp, arguments->ofi);
  gt_option_parser_set_mailaddress(op, "<steinbiss@zbh.uni-hamburg.de>");

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  gt_option_parser_refer_to_manual(op);
  gt_option_parser_set_min_max_args(op, 0, 1);

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
  GtGenomeNode *gn;
  GtRegionMapping *regionmapping;

  int had_err      = 0,
      tests_to_run = 0,
      arg = parsed_args;
  gt_error_check(err);
  gt_assert(arguments);

#ifdef HAVE_HMMER
  /* set additional arguments/options */
  arguments->pdom_opts.thresh.globT   = -FLT_MAX;
  arguments->pdom_opts.thresh.domT    = -FLT_MAX;
  arguments->pdom_opts.thresh.domE    = FLT_MAX;
  arguments->pdom_opts.thresh.autocut = CUT_NONE;
  arguments->pdom_opts.thresh.Z       = 1;
  arguments->pdom_opts.thresh.globE   = arguments->pdom_opts.evalue_cutoff;
#endif

  /* create region mapping */
  regionmapping = gt_seqid2file_regionmapping_new(arguments->seqfile,
                                                  arguments->regionmapping,
                                                  err);
  if (!regionmapping)
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
    arguments->pdom_opts.plan7_ts = gt_array_new(sizeof (struct plan7_s*));
    had_err = gt_pdom_load_hmm_files(&arguments->pdom_opts, err);
  }
#endif

  if (!had_err)
  {
    /* set up stream flow
     * ------------------*/
    last_stream = gff3_in_stream  = gt_gff3_in_stream_new_sorted(argv[arg]);

    last_stream = ltrdigest_stream = gt_ltrdigest_stream_new(last_stream,
                                                  tests_to_run,
                                                  regionmapping,
                                                  &arguments->pbs_opts,
                                                  &arguments->ppt_opts
#ifdef HAVE_HMMER
                                                 /*  */,&arguments->pdom_opts
#endif
                                      );

    /* attach tabular output stream, if requested */
    if (gt_str_length(arguments->prefix) > 0)
    {
      last_stream = tab_out_stream = gt_ltr_fileout_stream_new(last_stream,
                                              tests_to_run,
                                              regionmapping,
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
    }

    last_stream = gff3_out_stream = gt_gff3_out_stream_new(last_stream,
                                                           arguments->outfp);

    /* pull the features through the stream and free them afterwards */
    while (!(had_err = gt_node_stream_next(last_stream, &gn, err)) &&
           gn)
    {
      gt_genome_node_delete(gn);
    }
    gt_node_stream_delete(gff3_out_stream);
    gt_node_stream_delete(ltrdigest_stream);
    if (tab_out_stream)
      gt_node_stream_delete(tab_out_stream);
    gt_node_stream_delete(gff3_in_stream);
#ifdef HAVE_HMMER
    gt_pdom_clear_hmms(arguments->pdom_opts.plan7_ts);
#endif
  }
#ifdef HAVE_HMMER
  else
     gt_array_delete(arguments->pdom_opts.plan7_ts);
#endif

  gt_region_mapping_delete(regionmapping);
  gt_bioseq_delete(arguments->pbs_opts.trna_lib);

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
