/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/option.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "extended/gff3_out_stream.h"
#include "extended/bed_in_stream_api.h"
#include "tools/gt_bed_to_gff3.h"

static GtOptionParser* gt_bed_to_gff3_option_parser_new(GT_UNUSED
                                                        void *tool_arguments)
{
  GtOptionParser *op;
  op = gt_option_parser_new("[bed_file]",
                            "Parse BED file and show it as GFF3.");
  gt_option_parser_set_max_args(op, 1);
  return op;
}

static int gt_bed_to_gff3_runner(GT_UNUSED int argc, const char **argv,
                                 int parsed_args,
                                 GT_UNUSED void *tool_arguments, GtError *err)
{
  GtNodeStream *bed_in_stream = NULL, *gff3_out_stream = NULL;
  GtGenomeNode *gn;
  int had_err;

  gt_error_check(err);

  /* create a BED input stream */
  bed_in_stream = gt_bed_in_stream_new(argv[parsed_args]);

  /* create a GFF3 output stream */
  /* XXX: use proper genfile */
  gff3_out_stream = gt_gff3_out_stream_new(bed_in_stream, NULL);

  /* pull the features through the stream and free them afterwards */
  while (!(had_err = gt_node_stream_next(gff3_out_stream, &gn, err)) && gn)
    gt_genome_node_rec_delete(gn);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(bed_in_stream);

  return had_err;
}

GtTool* gt_bed_to_gff3(void)
{
  return gt_tool_new(NULL,
                     NULL,
                     gt_bed_to_gff3_option_parser_new,
                     NULL,
                     gt_bed_to_gff3_runner);
}
