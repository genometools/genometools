/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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
#include "libgtcore/bioseq.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/outputfile.h"
#include "libgtcore/safearith.h"
#include "libgtcore/unused.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gff3_out_stream.h"
#include "libgtltr/ltrdigest_def.h"
#include "libgtltr/ltrdigest_stream.h"
#include "libgtltr/ltrfileout_stream.h"
#include "tools/gt_ltrdigest.h"

typedef struct LTRdigestOptions {
  PBSOptions  pbs_opts;
  PPTOptions  ppt_opts;
  PdomOptions pdom_opts;
  Str *trna_lib;
  Str *prefix;
  bool verbose;
  OutputFileInfo *ofi;
  GenFile *outfp;
} LTRdigestOptions;

static void* gt_ltrdigest_arguments_new(void)
{
  LTRdigestOptions *arguments = ma_calloc(1, sizeof *arguments);
  memset(arguments, 0, sizeof *arguments);
  arguments->pdom_opts.hmm_files = strarray_new();
  arguments->trna_lib = str_new();
  arguments->prefix = str_new();
  arguments->ofi = outputfileinfo_new();
  return arguments;
}

static void gt_ltrdigest_arguments_delete(void *tool_arguments)
{
  LTRdigestOptions *arguments = tool_arguments;
  if (!arguments) return;
  strarray_delete(arguments->pdom_opts.hmm_files);
  str_delete(arguments->trna_lib);
  str_delete(arguments->prefix);
  genfile_close(arguments->outfp);
  outputfileinfo_delete(arguments->ofi);
  ma_free(arguments);
}

static OptionParser* gt_ltrdigest_option_parser_new(void *tool_arguments)
{
  LTRdigestOptions *arguments = tool_arguments;
  OptionParser *op;
  Option *o, *ot, *oh, *oto;
  static Range pptlen_defaults           = { 8, 30},
               uboxlen_defaults          = { 3, 30},
               pbsalilen_defaults        = {11, 30},
               pbsoffsetlen_defaults     = { 0,  5},
               pbstrnaoffsetlen_defaults = { 0,  5};
  assert(arguments);

  /* init */
  op = option_parser_new("[option ...] gff3_file sequence_file",
                         "Discovers and annotates sequence features in LTR "
                         "retrotransposon candidates.");

  /* PPT search options */

/*  o = option_new_uint("pptminlen",
                      "minimal required PPT length",
                      &arguments->ppt_opts.ppt_minlen,
                      6);
  option_parser_add_option(op, o); */

  o = option_new_range("pptlen",
                       "required PPT length range",
                       &arguments->ppt_opts.ppt_len,
                       &pptlen_defaults);
  option_parser_add_option(op, o);

  /*o = option_new_uint("uboxminlen",
                      "minimal required U-box length",
                      &arguments->ppt_opts.ubox_minlen,
                      3);
  option_parser_add_option(op, o); */

  o = option_new_range("uboxlen",
                       "required U-box length range",
                       &arguments->ppt_opts.ubox_len,
                       &uboxlen_defaults);
  option_parser_add_option(op, o);

  o = option_new_uint("pptradius",
                      "radius around beginning of 3' LTR "
                      "to search for PPT",
                      &arguments->ppt_opts.radius,
                      30);
  option_parser_add_option(op, o);

  /* PBS search options */

  ot = option_new_filename("trnas",
                          "tRNA library in multiple FASTA format for PBS "
                          "detection\n"
                          "Omit this option to disable PBS search.",
                          arguments->trna_lib);
  option_parser_add_option(op, ot);
  option_hide_default(ot);

  /* o = option_new_uint("pbsaliminlen",
                      "minimal required length of PBS/tRNA alignments",
                      &arguments->pbs_opts.ali_min_len,
                      11);
  option_parser_add_option(op, o);
  option_imply(o, ot); */

  o = option_new_range("pbsalilen",
                       "required PBS/tRNA alignment length range",
                       &arguments->pbs_opts.alilen,
                       &pbsalilen_defaults);
  option_parser_add_option(op, o);

  o = option_new_range("pbsoffset",
                       "allowed PBS offset from LTR boundary range",
                       &arguments->pbs_opts.offsetlen,
                       &pbsoffsetlen_defaults);
  option_parser_add_option(op, o);
  option_imply(o, ot);

  o = option_new_range("pbstrnaoffset",
                       "allowed PBS/tRNA 3' end alignment offset range",
                       &arguments->pbs_opts.trnaoffsetlen,
                       &pbstrnaoffsetlen_defaults);
  option_parser_add_option(op, o);
  option_imply(o, ot);

  o = option_new_uint("pbsmaxedist",
                      "maximal allowed PBS/tRNA alignment unit edit distance",
                      &arguments->pbs_opts.max_edist,
                      1);
  option_parser_add_option(op, o);
  option_imply(o, ot);

  o = option_new_uint("pbsradius",
                      "radius around end of 5' LTR "
                      " to search for PBS",
                      &arguments->pbs_opts.radius,
                      30);
  option_parser_add_option(op, o);
  option_imply(o, ot);

 /* Protein domain search options */

  oh = option_new_filenamearray("hmms",
                               "profile HMM models for domain detection "
                               "(separate by spaces, finish with --) in HMMER"
                               "2 format\n"
                               "Omit this option to disable pHMM search.",
                               arguments->pdom_opts.hmm_files);
  option_parser_add_option(op, oh);

  o = option_new_probability("pdomevalcutoff",
                             "E-value cutoff for pHMM search",
                             &arguments->pdom_opts.evalue_cutoff,
                             0.000001);
  option_parser_add_option(op, o);
  option_is_extended_option(o);
  option_imply(o, oh);

  o = option_new_uint_min("threads",
                          "number of concurrent worker threads to use in "
                          "pHMM scanning",
                          &arguments->pdom_opts.nof_threads,
                          2, 1);
  option_parser_add_option(op, o);

  o = option_new_uint("maxgaplen",
                      "maximal allowed gap size between fragments (in amino "
                      "acids) when chaining pHMM hits for a protein domain",
                      &arguments->pdom_opts.chain_max_gap_length,
                      50);
  option_parser_add_option(op, o);
  option_is_extended_option(o);

  /* Extended PBS options */

  o = option_new_int("pbsmatchscore",
                     "match score for PBS/tRNA alignments",
                      &arguments->pbs_opts.ali_score_match,
                      5);
  option_parser_add_option(op, o);
  option_is_extended_option(o);
  option_imply(o, ot);

  o = option_new_int("pbsmismatchscore",
                     "mismatch score for PBS/tRNA alignments",
                      &arguments->pbs_opts.ali_score_mismatch,
                      -10);
  option_parser_add_option(op, o);
  option_is_extended_option(o);
  option_imply(o, ot);

  o = option_new_int("pbsinsertionscore",
                     "insertion score for PBS/tRNA alignments",
                      &arguments->pbs_opts.ali_score_insertion,
                      -20);
  option_parser_add_option(op, o);
  option_is_extended_option(o);
  option_imply(o, ot);

  o = option_new_int("pbsdeletionscore",
                     "deletion score for PBS/tRNA alignments",
                      &arguments->pbs_opts.ali_score_deletion,
                      -20);
  option_parser_add_option(op, o);
  option_is_extended_option(o);
  option_imply(o, ot);

  /* Output files */
  oto = option_new_string("outfileprefix",
                          "prefix for output files (e.g. 'foo' will create "
                          "files called 'foo_*.csv' and 'foo_*.fas')\n"
                          "Omit this option for GFF3 output only.",
                          arguments->prefix,
                          NULL);
  option_parser_add_option(op, oto);
  option_hide_default(oto);

  /* verbosity */
  o = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, o);

  /* output file options */
  outputfile_register_options(op, &arguments->outfp, arguments->ofi);
  option_parser_set_mailaddress(op, "<ssteinbiss@stud.zbh.uni-hamburg.de>");

  option_parser_set_min_max_args(op, 2, 2);

  return op;
}

int gt_ltrdigest_arguments_check(UNUSED int rest_argc, void *tool_arguments,
                                 Error* err)
{
  LTRdigestOptions *arguments = tool_arguments;
  int had_err  = 0;

  /* TODO: more checks */
  /* -trnas */
  if (arguments->trna_lib
        && str_length(arguments->trna_lib) > 0)
  {
    if(!file_exists(str_get(arguments->trna_lib)))
    {
      error_set(err, "File '%s' does not exist!", str_get(arguments->trna_lib));
      had_err = -1;
    }
  }

  return had_err;
}

static int gt_ltrdigest_runner(UNUSED int argc, UNUSED const char **argv,
                               void *tool_arguments, Error *err)
{
  LTRdigestOptions *arguments = tool_arguments;
  GenomeStream *gff3_in_stream   = NULL,
               *gff3_out_stream  = NULL,
               *ltrdigest_stream = NULL,
               *tab_out_stream   = NULL,
               *last_stream      = NULL;
  GenomeNode *gn;

  int had_err      = 0,
      tests_to_run = 0;
  error_check(err);
  assert(arguments);

  /* set additional arguments/options */
  arguments->pdom_opts.thresh.globT   = -FLT_MAX;
  arguments->pdom_opts.thresh.domT    = -FLT_MAX;
  arguments->pdom_opts.thresh.domE    = FLT_MAX;
  arguments->pdom_opts.thresh.autocut = CUT_NONE;
  arguments->pdom_opts.thresh.Z       = 1;
  arguments->pdom_opts.thresh.globE   = arguments->pdom_opts.evalue_cutoff;

  /* Open sequence file */
  Bioseq *bioseq = bioseq_new(argv[1], err);
  if (error_is_set(err))
    had_err = -1;

  /* Always search for PPT. */
  tests_to_run |= LTRDIGEST_RUN_PPT;

  /* Open tRNA library if given. */
  if (!had_err && arguments->trna_lib
        && str_length(arguments->trna_lib) > 0)
  {
    tests_to_run |= LTRDIGEST_RUN_PBS;
    arguments->pbs_opts.trna_lib = bioseq_new(str_get(arguments->trna_lib),err);
    if (error_is_set(err))
      had_err = -1;
  }
  /* Open HMMER files if given. */
  if (!had_err && strarray_size(arguments->pdom_opts.hmm_files) > 0)
  {
    tests_to_run |= LTRDIGEST_RUN_PDOM;
    arguments->pdom_opts.plan7_ts = array_new(sizeof (struct plan7_s*));
    had_err = pdom_load_hmm_files(&arguments->pdom_opts,
                                  err);
  }

  if(!had_err)
  {
    /* set up stream flow
     * ------------------*/
    gff3_in_stream  = gff3_in_stream_new_sorted(argv[0],
                                                arguments->verbose &&
                                                arguments->outfp);

    ltrdigest_stream = ltrdigest_stream_new(gff3_in_stream,
                                            tests_to_run,
                                            bioseq,
                                            &arguments->pbs_opts,
                                            &arguments->ppt_opts,
                                            &arguments->pdom_opts);

    /* attach tabular output stream, if requested */
    if (str_length(arguments->prefix) > 0)
    {
      tab_out_stream = ltr_fileout_stream_new(ltrdigest_stream,
                                              tests_to_run,
                                              bioseq,
                                              str_get(arguments->prefix),
                                              &arguments->ppt_opts,
                                              &arguments->pbs_opts,
                                              &arguments->pdom_opts,
                                              str_get(arguments->trna_lib),
                                              argv[1],
                                              argv[0]);
      last_stream = tab_out_stream;
    }
    else
    {
      last_stream = ltrdigest_stream;
     }

    gff3_out_stream = gff3_out_stream_new(last_stream, arguments->outfp);

    /* pull the features through the stream and free them afterwards */
    while (!(had_err = genome_stream_next_tree(gff3_out_stream, &gn, err)) &&
           gn) {
      genome_node_rec_delete(gn);
    }

    genome_stream_delete(gff3_out_stream);
    genome_stream_delete(ltrdigest_stream);
    if (tab_out_stream)
    {
      genome_stream_delete(tab_out_stream);
    }
    genome_stream_delete(gff3_in_stream);
    pdom_clear_hmms(arguments->pdom_opts.plan7_ts);
  }

  bioseq_delete(bioseq);
  bioseq_delete(arguments->pbs_opts.trna_lib);

  return had_err;
}

Tool* gt_ltrdigest(void)
{
  return tool_new(gt_ltrdigest_arguments_new,
                  gt_ltrdigest_arguments_delete,
                  gt_ltrdigest_option_parser_new,
                  gt_ltrdigest_arguments_check,
                  gt_ltrdigest_runner);
}
