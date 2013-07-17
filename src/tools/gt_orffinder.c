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
#include "core/safearith.h"
#include "core/str_array.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/region_mapping.h"
#include "extended/orf_finder_stream.h"
#include "extended/seqid2file.h"
#include "gt_orffinder.h"

typedef struct {
  GtFile *outfp;
  unsigned int min,
               max;
  GtStrArray *types;
  GtOutputFileInfo *ofi;
  bool allorfs,
       verbose;
  GtSeqid2FileInfo *s2fi;
} GtOrffinderArguments;

static void* gt_orffinder_arguments_new(void)
{
  GtOrffinderArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  arguments->types = gt_str_array_new();
  arguments->s2fi = gt_seqid2file_info_new();
  return arguments;
}

static void gt_orffinder_arguments_delete(void *tool_arguments)
{
  GtOrffinderArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_str_array_delete(arguments->types);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_free(arguments);
}

static GtOptionParser* gt_orffinder_option_parser_new(void *tool_arguments)
{
  GtOrffinderArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optiontypes;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [indexname] [GFF3_file ...]",
                            "Identifies ORFs (open reading frames) in "
                            "sequences.");

  /* types */
  optiontypes = gt_option_new_string_array("types",
                                           "Specify regions which should be "
                                           "searched for open reading frames,\n"
                                           "e.g. 'LTR_retrotransposon'",
                                           arguments->types);
  gt_option_parser_add_option(op, optiontypes);

  /* -allorfs */
  option = gt_option_new_bool("allorfs", "search for all ORFs "
                              "instead of only the longest",
                              &arguments->allorfs, false);
  gt_option_parser_add_option(op, option);

  /* -min */
  option = gt_option_new_uint_min("min", "minimum length of ORF",
                                  &arguments->min, 30, 30);
  gt_option_parser_add_option(op, option);

  /* -max */
  option = gt_option_new_uint_max("max", "maximum length of ORF",
                                  &arguments->max, 10000, 1000000);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  /* region mapping and sequence source options */

  gt_seqid2file_register_options_ext(op, arguments->s2fi, false, false);

  return op;
}

static int gt_orffinder_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments, GtError *err)
{
  GtOrffinderArguments *arguments = tool_arguments;
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

static int gt_orffinder_runner(int argc, const char **argv,
                               int parsed_args, void *tool_arguments,
                               GtError *err)
{
  GtOrffinderArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream   = NULL,
               *last_stream      = NULL,
               *orffinder_stream = NULL,
               *gff3_out_stream  = NULL;
  int had_err = 0, i = 0, arg = parsed_args;
  GtLogger *logger = gt_logger_new(arguments->verbose,
                                   GT_LOGGER_DEFLT_PREFIX, stderr);
  GtHashmap *types = NULL;
  GtRegionMapping *rmap = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  types = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);

  if (gt_str_array_size(arguments->types) == 0) {
    gt_hashmap_add(types, (void*) "all", (void*) 1);
  } else {
    for (i = 0; i < gt_str_array_size(arguments->types); i++) {
      gt_hashmap_add(types,
                     (void*) gt_str_array_get(arguments->types, i),
                     (void*) 1);
      gt_hashmap_add(types, (void*) "all", (void*) 0);
    }
  }

  /* create region mapping */
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
      gt_error_set(err, "missing argument(s)");
      had_err = -1;
    }
    if (!had_err) {
      el = gt_encseq_loader_new();
      gt_encseq_loader_disable_autosupport(el);
      gt_encseq_loader_require_md5_support(el);
      gt_encseq_loader_require_description_support(el);
      encseq = gt_encseq_loader_load(el, argv[arg], err);
      arg++;
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

  if (!had_err) {
    last_stream = gff3_in_stream  = gt_gff3_in_stream_new_unsorted(argc - arg,
                                                                   argv + arg);
    last_stream = orffinder_stream = gt_orf_finder_stream_new(last_stream,
                                                             rmap,
                                                             types,
                                                             arguments->min,
                                                             arguments->max,
                                                             arguments->allorfs,
                                                             err);
    if (!orffinder_stream)
      had_err = -1;

    if (!had_err)  {
      last_stream = gff3_out_stream = gt_gff3_out_stream_new(last_stream,
                                                             arguments->outfp);

      /* pull the features through the stream and free them afterwards */
      had_err = gt_node_stream_pull(last_stream, err);
    }

    gt_node_stream_delete(orffinder_stream);
    gt_node_stream_delete(gff3_in_stream);
    gt_node_stream_delete(gff3_out_stream);
  }
  gt_hashmap_delete(types);
  types = NULL;
  gt_logger_delete(logger);
  gt_region_mapping_delete(rmap);

  return had_err;
}

GtTool* gt_orffinder(void)
{
  return gt_tool_new(gt_orffinder_arguments_new,
                     gt_orffinder_arguments_delete,
                     gt_orffinder_option_parser_new,
                     gt_orffinder_arguments_check,
                     gt_orffinder_runner);
}
