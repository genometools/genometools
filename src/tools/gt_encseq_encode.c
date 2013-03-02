/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include "core/alphabet.h"
#include "core/basename_api.h"
#include "core/ma.h"
#include "core/encseq.h"
#include "core/encseq_options.h"
#include "core/fileutils.h"
#include "core/logger_api.h"
#include "core/str_array_api.h"
#include "core/unused_api.h"
#include "tools/gt_encseq_encode.h"

typedef struct {
  GtEncseqOptions *eopts;
  bool showstats,
       no_esq_header,
       verbose;
  GtStr *indexname;
} GtEncseqEncodeArguments;

static void* gt_encseq_encode_arguments_new(void)
{
  GtEncseqEncodeArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->indexname = gt_str_new();
  return arguments;
}

static void gt_encseq_encode_arguments_delete(void *tool_arguments)
{
  GtEncseqEncodeArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_encseq_options_delete(arguments->eopts);
  gt_str_delete(arguments->indexname);
  gt_free(arguments);
}

static GtOptionParser* gt_encseq_encode_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option;
  GtEncseqEncodeArguments *arguments =
                                      (GtEncseqEncodeArguments*) tool_arguments;

  /* init */
  op = gt_option_parser_new("sequence_file [sequence_file "
                            "[sequence_file ...]]",
                            "Encode sequence files efficiently.");

  /* -showstats */
  option = gt_option_new_bool("showstats",
                              "show compression results",
                              &arguments->showstats,
                              false);
  gt_option_parser_add_option(op, option);

  /* -no_esq_header */
  option = gt_option_new_bool("no_esq_header",
                              "omit the header in the .esq file",
                              &arguments->no_esq_header,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* encoded sequence options */
  arguments->eopts = gt_encseq_options_register_encoding(op,
                                                         arguments->indexname,
                                                         NULL);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_args(op, 1);

  return op;
}

static int encode_sequence_files(GtStrArray *infiles, GtEncseqOptions *opts,
                                 const char *indexname, bool verbose,
                                 bool esq_no_header,
                                 GtError *err)
{
  GtEncseqEncoder *encseq_encoder;
  GtLogger *logger;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(infiles && gt_str_array_size(infiles) > 0 && opts);
  logger = gt_logger_new(verbose, "# ", stderr);
  encseq_encoder = gt_encseq_encoder_new_from_options(opts, err);
  if (!encseq_encoder)
    had_err = -1;
  if (!had_err) {
    gt_encseq_encoder_set_logger(encseq_encoder, logger);
    if (esq_no_header)
    {
      gt_encseq_encoder_disable_esq_header(encseq_encoder);
    }
    had_err = gt_encseq_encoder_encode(encseq_encoder, infiles, indexname, err);
  }
  gt_encseq_encoder_delete(encseq_encoder);
  gt_logger_delete(logger);
  return had_err;
}

static off_t index_size(const char *prefix, const char *suffix)
{
 off_t size = 0;
 GtStr *index_name = gt_str_new_cstr(prefix);
 gt_str_append_cstr(index_name, suffix);
 if (gt_file_exists(gt_str_get(index_name))) {
   size = gt_file_size(gt_str_get(index_name));
 }
 gt_str_delete(index_name);
 return size;
}

static void show_encoded_statistics(GtStrArray *infiles, const char *indexname)
{
  int i;
  off_t orig_size = 0, enc_size = 0;
  const char *seqfile;
  gt_assert(infiles);
  for (i=0; i < gt_str_array_size(infiles);i ++) {
    seqfile = gt_str_array_get(infiles, i);
    orig_size += gt_file_size(seqfile);
  }
  enc_size += index_size(indexname, GT_ALPHABETFILESUFFIX);
  enc_size += index_size(indexname, GT_ENCSEQFILESUFFIX);
  enc_size += index_size(indexname, GT_SSPTABFILESUFFIX);
  enc_size += index_size(indexname, GT_DESTABFILESUFFIX);
  enc_size += index_size(indexname, GT_SDSTABFILESUFFIX);
  enc_size += index_size(indexname, GT_OISTABFILESUFFIX);
  printf("encoded sequence file(s) are %.1f%% of original file size\n",
         ((double) enc_size / orig_size) * 100.0);
}

static int gt_encseq_encode_runner(GT_UNUSED int argc, const char **argv,
                               int parsed_args, GT_UNUSED void *tool_arguments,
                               GtError *err)
{
  int had_err = 0,
      i;
  GtEncseqEncodeArguments *arguments =
                                      (GtEncseqEncodeArguments*) tool_arguments;
  GtStrArray *infiles;
  gt_error_check(err);

  infiles = gt_str_array_new();
  for (i = parsed_args; i < argc; i++) {
    gt_str_array_add_cstr(infiles, argv[i]);
  }

  if (gt_str_length(arguments->indexname) == 0UL) {
    if (gt_str_array_size(infiles) > 1UL) {
      gt_error_set(err,"if more than one input file is given, then "
                       "option -indexname is mandatory");
      had_err = -1;
    } else {
      char *basenameptr;
      basenameptr = gt_basename(gt_str_array_get(infiles, 0UL));
      gt_str_set(arguments->indexname, basenameptr);
      gt_free(basenameptr);
    }
  }

  if (!had_err) {
    gt_assert(gt_str_length(arguments->indexname) > 0UL);
    had_err = encode_sequence_files(infiles,
                                    arguments->eopts,
                                    gt_str_get(arguments->indexname),
                                    arguments->verbose,
                                    arguments->no_esq_header,
                                    err);
  }

  if (!had_err && arguments->showstats)
    show_encoded_statistics(infiles, gt_str_get(arguments->indexname));

  gt_str_array_delete(infiles);
  return had_err;
}

GtTool* gt_encseq_encode(void)
{
  return gt_tool_new(gt_encseq_encode_arguments_new,
                     gt_encseq_encode_arguments_delete,
                     gt_encseq_encode_option_parser_new,
                     NULL,
                     gt_encseq_encode_runner);
}
