  /*
  Copyright (c) 2010      Sascha Kastens <mail@skastens.de>
  Copyright (c) 2011-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2013 Center for Bioinformatics, University of Hamburg

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
#include "core/codon_iterator_encseq_api.h"
#include "core/encseq.h"
#include "core/fileutils_api.h"
#include "core/hashmap.h"
#include "core/log.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/str_array.h"
#include "core/output_file_api.h"
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/str_api.h"
#include "core/str_array.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/region_mapping.h"
#include "extended/orf_finder_stream.h"
#include "extended/seqid2file.h"
#include "gt_orfgenerator.h"

typedef struct {
  GtFile *outfp;
  unsigned int min,
               max;
  GtStr *type;
  bool all;
  GtOutputFileInfo *ofi;
} GtorfgeneratorArguments;

static void* gt_orfgenerator_arguments_new(void)
{
  GtorfgeneratorArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  arguments->type = gt_str_new();
  return arguments;
}

static void gt_orfgenerator_arguments_delete(void *tool_arguments)
{
  GtorfgeneratorArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_str_delete(arguments->type);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_orfgenerator_option_parser_new(void *tool_arguments) {
  GtorfgeneratorArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [sequence_file ...]",
                            "Identifies ORFs (open reading frames) in "
                            "sequences.");

  /* -min */
  option = gt_option_new_uint_min("min", "minimum length of ORF",
                                  &arguments->min, 30, 30);
  gt_option_parser_add_option(op, option);

  /* -max */
  option = gt_option_new_uint_max("max", "maximum length of ORF",
                                  &arguments->max, 10000, 1000000);
  gt_option_parser_add_option(op, option);

  /* -type */
  option = gt_option_new_string("type", "type for ORF output",
                                arguments->type, "CDS");
  gt_option_parser_add_option(op, option);

  /* -all */
  option = gt_option_new_bool("all", "output all overlapping ORFs",
                              &arguments->all, false);
  gt_option_parser_add_option(op, option);

  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  gt_option_parser_set_min_args(op, 1);

  return op;
}

static int gt_orfgenerator_arguments_check(GT_UNUSED int rest_argc,
                                           void *tool_arguments, GtError *err)
{
  GtorfgeneratorArguments *arguments = tool_arguments;
  int had_err = 0;
  int diff;

  gt_error_check(err);
  gt_assert(arguments);

  diff = (arguments->max - arguments->min);
  if (diff < 0) {
    gt_error_set(err, "Value for -min must be larger than -max");
    had_err = -1;
  }

  return had_err;
}

static int gt_orfgenerator_runner(int argc, const char **argv,
                               int parsed_args, void *tool_arguments,
                               GtError *err)
{
  GtorfgeneratorArguments *arguments = tool_arguments;
  GtStrArray *filenames;
  GtSeqIterator *seqit;
  GtNodeStream *orfgenerator_stream = NULL,
               *gff3_out_stream = NULL;
  int had_err = 0, i = 0, arg = parsed_args;

  gt_error_check(err);
  gt_assert(arguments);

  filenames = gt_str_array_new();

  for (i = arg; i < argc; i++) {
    gt_str_array_add_cstr(filenames, argv[i]);
  }

  seqit = gt_seq_iterator_sequence_buffer_new(filenames, err);
  if (!seqit)
    had_err = -1;

  if (!had_err) {
    orfgenerator_stream = gt_orf_finder_stream_new_from_seq(seqit,
                                                            arguments->min,
                                                            arguments->max,
                                                            arguments->all,
                                                            err);
    if (!orfgenerator_stream)
      had_err = -1;
  }

  if (!had_err) {
    gt_orf_finder_stream_set_type((GtORFFinderStream*) orfgenerator_stream,
                                  gt_str_get(arguments->type));
    gff3_out_stream = gt_gff3_out_stream_new(orfgenerator_stream, arguments->outfp);
  }

  if (!had_err)
    had_err = gt_node_stream_pull(gff3_out_stream, err);

  gt_seq_iterator_delete(seqit);
  gt_str_array_delete(filenames);
  gt_node_stream_delete(orfgenerator_stream);
  gt_node_stream_delete(gff3_out_stream);

  return had_err;
}

GtTool* gt_orfgenerator(void)
{
  return gt_tool_new(gt_orfgenerator_arguments_new,
                     gt_orfgenerator_arguments_delete,
                     gt_orfgenerator_option_parser_new,
                     gt_orfgenerator_arguments_check,
                     gt_orfgenerator_runner);
}
