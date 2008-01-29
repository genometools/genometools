/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/option.h"
#include "libgtcore/outputfile.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/chseqids_stream.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gff3_out_stream.h"
#include "libgtext/gtdatahelp.h"
#include "libgtext/sort_stream.h"
#include "tools/gt_chseqids.h"

#define DEFAULT_JOINLENGTH 300

typedef struct {
  bool sort,
       verbose;
  GenFile *outfp;
} ChseqidsArguments;

static OPrval parse_options(int *parsed_args, ChseqidsArguments *arguments,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  OutputFileInfo *ofi;
  Option *option;
  OPrval oprval;
  error_check(err);

  /* init */
  op = option_parser_new("[option ...] mapping_file [GFF3_file]",
                         "Change sequence ids by the mapping given in "
                         "mapping_file.");
  ofi = outputfileinfo_new();

  /* -sort */
  option = option_new_bool("sort", "sort the GFF3 features after changing the "
                           "sequence ids\n(memory consumption is O(file_size))",
                           &arguments->sort, false);
  option_parser_add_option(op, option);

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* output file options */
  outputfile_register_options(op, &arguments->outfp, ofi);

  /* parse options */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  option_parser_set_min_max_args(op, 1, 2);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);

  /* free */
  outputfileinfo_delete(ofi);
  option_parser_delete(op);

  return oprval;
}

int gt_chseqids(int argc, const char **argv, Error *err)
{
  GenomeStream *gff3_in_stream, *chseqids_stream, *sort_stream = NULL,
               *gff3_out_stream = NULL;
  GenomeNode *gn;
  ChseqidsArguments arguments;
  Str *chseqids;
  int parsed_args, had_err = 0;

  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create the streams */
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args + 1],
                                             arguments.verbose &&
                                             arguments.outfp);
  chseqids = str_new_cstr(argv[parsed_args]);
  if (!(chseqids_stream = chseqids_stream_new(gff3_in_stream, chseqids, err)))
    had_err = -1;
  str_delete(chseqids);
  if (!had_err) {
    if (arguments.sort) {
      sort_stream = sort_stream_new(chseqids_stream);
      gff3_out_stream = gff3_out_stream_new(sort_stream, arguments.outfp);
    }
    else
      gff3_out_stream = gff3_out_stream_new(chseqids_stream, arguments.outfp);
  }

  /* pull the features through the stream and free them afterwards */
  if (!had_err) {
    while (!(had_err = genome_stream_next_tree(gff3_out_stream, &gn, err)) &&
           gn) {
      genome_node_rec_delete(gn);
    }
  }

  /* free */
  genome_stream_delete(gff3_out_stream);
  genome_stream_delete(chseqids_stream);
  genome_stream_delete(sort_stream);
  genome_stream_delete(gff3_in_stream);
  genfile_close(arguments.outfp);

  return had_err;
}
