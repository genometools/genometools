/*
  Copyright (c) 2010-2011 Gordon Gremme <gordon@gremme.org>

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

#include "core/compat.h"
#include "core/cstr_array.h"
#include "core/ma_api.h"
#include "core/output_file.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/yarandom.h"
#include "core/warning_api.h"
#include "extended/add_introns_stream_api.h"
#include "extended/buffer_stream.h"
#include "extended/feature_type.h"
#include "extended/gff3_in_stream_api.h"
#include "extended/seqid2file.h"
#include "extended/splice_site_info_stream.h"
#include "gth/bssm_train_stream.h"
#include "gth/gt_gthbssmtrain.h"

typedef struct {
  bool gcdonor,
       force,
       intermediate,
       verbose,
       gzip,
       bzip2;
  unsigned int seed;
  GtStr *outdir,
        *filter_type,
        *extract_type;
  unsigned int good_exon_count;
  double cutoff;
  GtFileMode file_mode;
  GtSeqid2FileInfo *s2fi;
} GthBSSMTrainArguments;

static void* gt_gthbssmtrain_arguments_new(void)
{
  GthBSSMTrainArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->outdir = gt_str_new();
  arguments->filter_type = gt_str_new();
  arguments->extract_type = gt_str_new();
  arguments->file_mode = GT_FILE_MODE_UNCOMPRESSED;
  arguments->s2fi = gt_seqid2file_info_new();
  return arguments;
}

static void gt_gthbssmtrain_arguments_delete(void *tool_arguments)
{
  GthBSSMTrainArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_str_delete(arguments->extract_type);
  gt_str_delete(arguments->filter_type);
  gt_str_delete(arguments->outdir);
  gt_free(arguments);
}

static GtOptionParser* gt_gthbssmtrain_option_parser_new(void *tool_arguments)
{
  GthBSSMTrainArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optgzip, *optbzip2;
  gt_assert(arguments);
  op = gt_option_parser_new("[option ...] GFF3_file", "Create BSSM training "
                            "data from annotation given in GFF3_file.");

  /* -outdir */
  option = gt_option_new_string("outdir", "set name of output directory to "
                                "which the training files are written",
                                arguments->outdir, "training_data");
  gt_option_parser_add_option(op, option);

  /* -gcdonor */
  option = gt_option_new_bool("gcdonor", "extract training data for GC donor "
                              "sites", &arguments->gcdonor, true);
  gt_option_parser_add_option(op, option);

  /* -filtertype */
  option = gt_option_new_string("filtertype", "set type of features to used "
                                "for filtering (usually 'exon' or 'CDS')",
                                arguments->filter_type, gt_ft_exon);
  gt_option_parser_add_option(op, option);

  /* -goodexoncount */
  option = gt_option_new_uint("goodexoncount", "set the minimum number of good "
                              "exons a feature must have to be included into "
                              "the training data", &arguments->good_exon_count,
                              1);
  gt_option_parser_add_option(op, option);

  /* -cutoff */
  option = gt_option_new_double("cutoff", "set the minimum score an exon must "
                                "have to count towards the ``good exon count'' "
                                "(exons without a score count as good)",
                                &arguments->cutoff, 1.0);
  gt_option_parser_add_option(op, option);

  /* -extracttype */
  option = gt_option_new_string("extracttype", "set type of features to be "
                                "extracted as exons (usually 'exon' or 'CDS')",
                                arguments->extract_type, gt_ft_CDS);
  gt_option_parser_add_option(op, option);

  /* -intermediate */
  option = gt_option_new_bool("intermediate", "write out files containing "
                              "intermediate results", &arguments->intermediate,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -seqfile, -matchdesc, -usedesc and -region_mapping */
  gt_seqid2file_register_options(op, arguments->s2fi);

  /* -seed */
  option = gt_option_new_uint("seed", "set seed for random number generator "
                              "manually\n0 generates a seed from the current "
                              "time and the process id", &arguments->seed, 0);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* -gzip */
  optgzip = gt_option_new_bool("gzip", "write gzip compressed output files",
                               &arguments->gzip, false);
  gt_option_parser_add_option(op, optgzip);

  /* -bzip2 */
  optbzip2 = gt_option_new_bool("bzip2", "write bzip2 compressed output files",
                                &arguments->bzip2, false);
  gt_option_parser_add_option(op, optbzip2);

  /* -force */
  option = gt_option_new_bool(GT_FORCE_OPT_CSTR, "force writing to output "
                              "files", &arguments->force, false);
  gt_option_parser_add_option(op, option);

  gt_option_exclude(optgzip, optbzip2);

  gt_option_parser_set_mail_address(op, "<gordon@gremme.org>");
  gt_option_parser_set_min_max_args(op, 1, 1);
  return op;
}

static int gt_gthbssmtrain_arguments_check(GT_UNUSED int rest_argc,
                                           void *tool_arguments,
                                           GT_UNUSED GtError *err)
{
  GthBSSMTrainArguments *arguments = tool_arguments;
  gt_error_check(err);
  gt_assert(arguments);
  if (arguments->gzip)
    arguments->file_mode = GT_FILE_MODE_GZIP;
  if (arguments->bzip2)
    arguments->file_mode = GT_FILE_MODE_BZIP2;
  return 0;
}

static GtBufferStream* get_buffer_and_show_splice_sites(const char *filename,
                                                        GtRegionMapping *rm,
                                                        bool show_gc,
                                                        GtError *err)
{
  GtBufferStream *buffer_stream_a, *buffer_stream_b;
  GtNodeStream *gff3_in_stream, *ssi_stream;
  bool processed_intron;
  int had_err;

  gt_error_check(err);

  gff3_in_stream = gt_gff3_in_stream_new_sorted(filename);

  buffer_stream_a = (GtBufferStream*) gt_buffer_stream_new(gff3_in_stream);
  ssi_stream = gt_splice_site_info_stream_new((GtNodeStream*) buffer_stream_a,
                                              gt_region_mapping_ref(rm));

  had_err = gt_node_stream_pull(ssi_stream, err);

  if (had_err) {
    gt_node_stream_delete(ssi_stream);
    gt_node_stream_delete((GtNodeStream*) buffer_stream_a);
    gt_node_stream_delete(gff3_in_stream);
    return NULL;
  }

  processed_intron = gt_splice_site_info_stream_intron_processed(ssi_stream);

  if (processed_intron) {
    if (!gt_splice_site_info_stream_show_canonical(ssi_stream, show_gc))
      gt_warning("no gt-ag or gc-ag splice sites found\n");
  }
  else {
    GtNodeStream *add_introns_stream;
    gt_node_stream_delete(ssi_stream);
    gt_buffer_stream_dequeue(buffer_stream_a);
    buffer_stream_b = (GtBufferStream*)
                      gt_buffer_stream_new((GtNodeStream*) buffer_stream_a);
    add_introns_stream = gt_add_introns_stream_new((GtNodeStream*)
                                                   buffer_stream_b);
    ssi_stream = gt_splice_site_info_stream_new(add_introns_stream,
                                                gt_region_mapping_ref(rm));
    had_err = gt_node_stream_pull(ssi_stream, err);
    if (had_err)
      gt_node_stream_delete((GtNodeStream*) buffer_stream_b);
    else {
      if (!gt_splice_site_info_stream_show_canonical(ssi_stream, show_gc))
        gt_warning("no gt-ag or gc-ag splice sites found\n");
    }
    gt_node_stream_delete(add_introns_stream);
    gt_node_stream_delete((GtNodeStream*) buffer_stream_a);
  }

  gt_node_stream_delete(ssi_stream);
  gt_node_stream_delete(gff3_in_stream);

  if (had_err)
    return NULL;
  return processed_intron ? buffer_stream_a : buffer_stream_b;
}

static int gt_gthbssmtrain_runner(GT_UNUSED int argc, const char **argv,
                                  int parsed_args, void *tool_arguments,
                                  GtError *err)
{
  GthBSSMTrainArguments *arguments = tool_arguments;
  GtNodeStream *bssm_train_stream = NULL;
  GtBufferStream *buffer_stream = NULL;
  GtRegionMapping *region_mapping;
  GthBSSMSeqProcessor *bsp = NULL;
  GtFile *logfp = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* create region mapping */
  region_mapping = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
  if (!region_mapping)
    had_err = -1;

  if (!had_err) {
    if (!(bsp = gth_bssm_seq_processor_new(gt_str_get(arguments->outdir),
                                           arguments->file_mode,
                                           arguments->force, arguments->gcdonor,
                                           err))) {
      had_err = -1;
    }
  }

  /* open log file */
  if (!had_err) {
    GtStr *logfile = gt_str_clone(arguments->outdir);
    gt_str_append_char(logfile, GT_PATH_SEPARATOR);
    gt_str_append_cstr(logfile, "gthbssmtrain.run");
    if (!(logfp = gt_output_file_xopen_forcecheck(gt_str_get(logfile), "w",
                                                  arguments->force,
                                                  err))) {
      had_err = -1;
    }
    gt_str_delete(logfile);
  }

  if (!had_err) {
    /* log argv */
    gt_file_xprintf(logfp, "arguments=");
    gt_cstr_array_show_genfile(argv + 1, logfp);

    /* initialize random number generator */
    arguments->seed = gt_ya_rand_init(arguments->seed);
    if (arguments->verbose)
      printf("seed=%u\n", arguments->seed);
    gt_file_xprintf(logfp, "seed=%u\n", arguments->seed);
  }

  if (!had_err) {
    buffer_stream = get_buffer_and_show_splice_sites(argv[parsed_args],
                                                     region_mapping,
                                                     arguments->gcdonor, err);
    if (!buffer_stream)
      had_err = -1;
    else
      gt_buffer_stream_dequeue(buffer_stream);
  }

  if (!had_err) {
    /* create BSSM training stream */
    bssm_train_stream =
      gth_bssm_train_stream_new((GtNodeStream*) buffer_stream, region_mapping,
                                bsp, gt_str_get(arguments->filter_type),
                                gt_str_get(arguments->extract_type),
                                arguments->good_exon_count, arguments->cutoff);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(bssm_train_stream, err);
  }

  if (!had_err) {
    gth_bssm_seq_processor_squash(bsp);
    had_err = gth_bssm_seq_processor_find_true_sites(bsp, region_mapping, err);
  }

  if (!had_err)
    had_err = gth_bssm_seq_processor_find_false_sites(bsp, region_mapping, err);

  if (!had_err && arguments->intermediate)
    had_err = gth_bssm_seq_processor_write_intermediate(bsp, err);

  if (!had_err) {
    gth_bssm_seq_processor_sample(bsp, arguments->verbose, logfp);
    gth_bssm_seq_processor_write(bsp);
  }

  gt_node_stream_delete(bssm_train_stream);
  gt_node_stream_delete((GtNodeStream*) buffer_stream);
  gt_file_delete(logfp);

  return had_err;
}

GtTool* gt_gthbssmtrain(void)
{
  return gt_tool_new(gt_gthbssmtrain_arguments_new,
                     gt_gthbssmtrain_arguments_delete,
                     gt_gthbssmtrain_option_parser_new,
                     gt_gthbssmtrain_arguments_check,
                     gt_gthbssmtrain_runner);
}
