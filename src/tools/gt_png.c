/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include <stdlib.h>
#include "cairo_stream.h"
#include "error.h"
#include "fileutils.h"
#include "gff3_in_stream.h"
#include "gff3_out_stream.h"
#include "option.h"
#include "showgthversionfunc.h"
#include "undef.h"

typedef struct {
  Str *sequence_region_id;
  unsigned long from,
                to;
  int width;
  unsigned int pipe,
               verbose;
} Png_info;

static int parse_options(Png_info *info, int argc, char **argv)
{
  int parsed_args;
  Option_parser *op = option_parser_new("[option ...] png_file [gff3_file]",
                                        "Draw gff3 features as png.");
  Option *option;
  unsigned int force;

  /* -sequence-region */
  option = option_new_string("sequence-region", "set drawed sequence-region "
                             "(by default the first one is taken)",
                             info->sequence_region_id, NULL);
  option_hide_default(option);
  option_parser_add_option(op, option);

  /* -from */
  option = option_new_ulong_min("from", "set position from where on features "
                                "are shown", &info->from, 1, 1);
  option_hide_default(option);
  option_parser_add_option(op, option);

  /* -to */
  option = option_new_ulong("to", "set maximal position of features to be "
                            "shown", &info->to, ~0UL);
  option_hide_default(option);
  option_parser_add_option(op, option);

  /* -width */
  option = option_new_int_min("width", "set the width of the png file",
                              &info->width, 1024, 1);
  option_parser_add_option(op, option);

  /* -pipe */
  option = option_new_boolean("pipe", "use pipe mode (i.e., show all gff3 "
                              "features on stdout)", &info->pipe, 0);
  option_parser_add_option(op, option);

  /* -force */
  option = option_new_boolean("force", "force writing to output file", &force,
                              0);
  option_parser_add_option(op, option);

  /* -v */
  option = option_new_verbose(&info->verbose);
  option_parser_add_option(op, option);

  /* parsing */
  parsed_args = option_parser_parse_min_max_args(op, argc, argv,
                                                 showgthversionfunc, 1, 2);

  /* checks */
  if (info->from > info->to)
    error("\"-from\" must be <= \"-to\"\n");

  if (!force && file_exists(argv[parsed_args]))
    error("file \"%s\" exists already. use option -force to overwrite",
          argv[parsed_args]);

  option_parser_free(op);
  return parsed_args;
}

int main(int argc, char *argv[])
{
  Genome_stream *gff3_in_stream,
                *cairo_stream,
                *gff3_out_stream = NULL;
  Genome_node *gn;
  int parsed_args;
  Png_info info;

  /* option parsing */
  info.sequence_region_id = str_new();
  parsed_args = parse_options(&info, argc, argv);

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args + 1],
                                             info.verbose && !info.pipe);

  /* create a cairo stream */
  cairo_stream = cairo_stream_new(gff3_in_stream, info.sequence_region_id,
                                  info.from, info.to, argv[parsed_args],
                                  info.width);

  /* create gff3 output stream if -pipe was used */
  if (info.pipe)
    gff3_out_stream = gff3_out_stream_new(cairo_stream, stdout);

  /* pull the features through the stream and free them afterwards */
  while ((gn = genome_stream_next_tree(info.pipe ? gff3_out_stream
                                                 : cairo_stream,
                                       NULL))) {
    genome_node_rec_free(gn);
  }

  /* draw */
  cairo_stream_draw((Cairo_stream*) cairo_stream, info.verbose, NULL);

  /* free */
  genome_stream_free(gff3_out_stream);
  genome_stream_free(cairo_stream);
  genome_stream_free(gff3_in_stream);
  str_free(info.sequence_region_id);

  return EXIT_SUCCESS;
}
