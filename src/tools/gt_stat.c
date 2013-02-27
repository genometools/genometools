/*
  Copyright (c) 2005-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "extended/add_introns_stream_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/sort_stream_api.h"
#include "extended/stat_stream_api.h"
#include "tools/gt_stat.h"

typedef struct {
  bool gene_length_distribution,
       gene_score_distribution,
       exon_number_distribution,
       exon_length_distribution,
       intron_length_distribution,
       cds_length_distribution,
       used_sources,
       addintrons,
       verbose;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} StatArguments;

static void* gt_stat_argument_new(void)
{
  StatArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_stat_arguments_delete(void *tool_arguments)
{
  StatArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_stat_option_parser_new(void *tool_arguments)
{
  StatArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;

  op = gt_option_parser_new("[option ...] [GFF3_file ...]",
                            "Show statistics about features contained in GFF3 "
                            "files.");

  /* -genelengthdistri */
  option = gt_option_new_bool("genelengthdistri",
                              "show gene length distribution",
                              &arguments->gene_length_distribution, false);
  gt_option_parser_add_option(op, option);

  /* -genescoresdistri */
  option = gt_option_new_bool("genescoredistri", "show gene score distribution",
                              &arguments->gene_score_distribution, false);
  gt_option_parser_add_option(op, option);

  /* -exonlengthdistri */
  option = gt_option_new_bool("exonlengthdistri",
                              "show exon length distribution",
                              &arguments->exon_length_distribution, false);
  gt_option_parser_add_option(op, option);

  /* -exonnumberdistri */
  option = gt_option_new_bool("exonnumberdistri",
                              "show exon number distribution",
                              &arguments->exon_number_distribution, false);
  gt_option_parser_add_option(op, option);

  /* -intronlengthdistri */
  option = gt_option_new_bool("intronlengthdistri",
                              "show intron length distribution",
                              &arguments->intron_length_distribution, false);
  gt_option_parser_add_option(op, option);

  /* -cdslengthdistri */
  option = gt_option_new_bool("cdslengthdistri", "show CDS length distribution "
                              "(measured in amino acids)",
                              &arguments->cds_length_distribution, false);
  gt_option_parser_add_option(op, option);

  /* -source */
  option = gt_option_new_bool("source", "show the set of used source tags "
                              "(column 2 in regular GFF3 lines)",
                              &arguments->used_sources, false);
  gt_option_parser_add_option(op, option);

  /* -addintrons */
  option = gt_option_new_bool("addintrons", "add intron features between "
                              "existing exon features (before computing stats)",
                              &arguments->addintrons, false);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  return op;
}

static int gt_stat_runner(int argc, const char **argv, int parsed_args,
                          void *tool_arguments, GtError *err)
{
  StatArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream, *sort_stream, *add_introns_stream = NULL,
               *stat_stream;
  int had_err;
  gt_error_check(err);

  /* create a gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);
  if (arguments->verbose)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);

  /* create add introns stream if -addintrons was used */
  if (arguments->addintrons) {
    sort_stream = gt_sort_stream_new(gff3_in_stream);
    add_introns_stream = gt_add_introns_stream_new(sort_stream);
  }

  /* create s status stream */
  stat_stream = gt_stat_stream_new(arguments->addintrons
                                   ? add_introns_stream : gff3_in_stream,
                                   arguments->gene_length_distribution,
                                   arguments->gene_score_distribution,
                                   arguments->exon_length_distribution,
                                   arguments->exon_number_distribution,
                                   arguments->intron_length_distribution,
                                   arguments->cds_length_distribution,
                                   arguments->used_sources);

  /* pull the features through the stream , compute the statistics, and free
     them afterwards */
  had_err = gt_node_stream_pull(stat_stream, err);

  /* show statistics */
  if (!had_err)
    gt_stat_stream_show_stats((GtStatStream*) stat_stream, arguments->outfp);

  /* free */
  gt_node_stream_delete(stat_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool *gt_stat(void)
{
  return gt_tool_new(gt_stat_argument_new,
                     gt_stat_arguments_delete,
                     gt_stat_option_parser_new,
                     NULL,
                     gt_stat_runner);
}
