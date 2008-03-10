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
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/outputfile.h"
#include "libgtcore/safearith.h"
#include "libgtcore/unused.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gff3_out_stream.h"
#include "libgtltr/pbs.h"
#include "libgtltr/ppt.h"
#include "libgtltr/pdom.h"
#include "tools/gt_ltrdigest.h"

typedef struct LTRdigestOptions {
  PBSOptions  pbs_opts;
  PPTOptions  ppt_opts;
  PdomOptions pdom_opts;
  Str *trna_lib;
  bool verbose;
  OutputFileInfo *ofi;
  GenFile *outfp;
} LTRdigestOptions;

static void* gt_ltrdigest_arguments_new(void)
{
  LTRdigestOptions *arguments = ma_calloc(1, sizeof *arguments);
  arguments->pdom_opts.hmm_files = strarray_new();
  arguments->trna_lib = str_new();
  arguments->ofi = outputfileinfo_new();
  return arguments;
}

static void gt_ltrdigest_arguments_delete(void *tool_arguments)
{
  LTRdigestOptions *arguments = tool_arguments;
  if (!arguments) return;
  strarray_delete(arguments->pdom_opts.hmm_files);
  str_delete(arguments->trna_lib);
  genfile_close(arguments->outfp);
  outputfileinfo_delete(arguments->ofi);
  ma_free(arguments);
}

static OptionParser* gt_ltrdigest_option_parser_new(void *tool_arguments)
{
  LTRdigestOptions *arguments = tool_arguments;
  OptionParser *op;
  Option *o;
  assert(arguments);

  /* init */
  op = option_parser_new("[option ...] gff3_file sequence_file",
                         "Discovers and annotates sequence features in LTR "
                         "retrotransposon candidates.");

  /* PPT search options */

  o = option_new_uint("pptminlen",
                      "minimum PPT length",
                      &arguments->ppt_opts.ppt_minlen,
                      6);
  option_parser_add_option(op, o);

  o = option_new_uint("uboxminlen",
                      "minimum U-box length",
                      &arguments->ppt_opts.ubox_minlen,
                      3);
  option_parser_add_option(op, o);

  o = option_new_uint("pptradius",
                      "radius around beginning of 3' LTR "
                      "to search for PPT",
                      &arguments->ppt_opts.radius,
                      30);
  option_parser_add_option(op, o);

  /* PBS search options */

  o = option_new_uint("pbsaliminlen",
                      "minimum length of PBS/tRNA alignments",
                      &arguments->pbs_opts.ali_min_len,
                      11);
  option_parser_add_option(op, o);

  o = option_new_uint("pbsmaxoffsetltr",
                      "maximal allowed PBS offset from LTR boundary",
                      &arguments->pbs_opts.max_offset,
                      5);
  option_parser_add_option(op, o);

  o = option_new_uint("pbsmaxoffsettrna",
                      "maximal allowed PBS alignment offset from tRNA 3' end",
                      &arguments->pbs_opts.max_offset_trna,
                      10);
  option_parser_add_option(op, o);

  o = option_new_uint("pbsmaxedist",
                      "maximal allowed PBS/tRNA alignment unit edit distance",
                      &arguments->pbs_opts.max_edist,
                      1);
  option_parser_add_option(op, o);

  o = option_new_uint("pbsradius",
                      "radius around end of 5' LTR "
                      " to search for PBS",
                      &arguments->pbs_opts.radius,
                      30);
  option_parser_add_option(op, o);

  o = option_new_filename("trnas",
                          "tRNA library in multiple FASTA format for PBS "
                          "detection",
                          arguments->trna_lib);
  option_parser_add_option(op, o);

 /* PBS search options */

  o = option_new_probability("pdomevalcutoff",
                             "E-value cutoff for HMMER protein domain search",
                             &arguments->pdom_opts.evalue_cutoff,
                             0.00001);
  option_parser_add_option(op, o);

  o = option_new_filenamearray("hmms",
                               "HMMER models for protein domain detection "
                               "(separate by spaces, finish with --) in HMMER "
                               "2.0 format",
                               arguments->pdom_opts.hmm_files);
  option_parser_add_option(op, o);

  /* verbosity */

  o = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, o);

  /* output file options */
  outputfile_register_options(op, &arguments->outfp, arguments->ofi);
  option_parser_set_mailaddress(op, "<ssteinbiss@stud.zbh.uni-hamburg.de>");

  option_parser_set_min_max_args(op, 2, 2);

  return op;
}

static int gt_ltrdigest_runner(UNUSED int argc, UNUSED const char **argv,
                               void *tool_arguments, Error *err)
{
  LTRdigestOptions *arguments = tool_arguments;
  GenomeStream *gff3_in_stream,
               *gff3_out_stream;
  GenomeNode *gn;

  int had_err = 0;
  error_check(err);
  assert(arguments);

  gff3_in_stream  = gff3_in_stream_new_sorted(argv[0],
                                              arguments->verbose &&
                                              arguments->outfp);

  /* insert LTRdigest stream here! */

  gff3_out_stream = gff3_out_stream_new(gff3_in_stream, arguments->outfp);

  /* pull the features through the stream and free them afterwards */
  while (!(had_err = genome_stream_next_tree(gff3_out_stream, &gn, err)) &&
         gn) {
    genome_node_rec_delete(gn);
  }

  genome_stream_delete(gff3_out_stream);
  genome_stream_delete(gff3_in_stream);
  return had_err;
}

Tool* gt_ltrdigest(void)
{
  return tool_new(gt_ltrdigest_arguments_new,
                  gt_ltrdigest_arguments_delete,
                  gt_ltrdigest_option_parser_new,
                  NULL,
                  gt_ltrdigest_runner);
}
