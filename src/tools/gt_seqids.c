/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/option_api.h"
#include "core/unused_api.h"
#include "extended/collect_ids_visitor.h"
#include "extended/gff3_in_stream.h"
#include "extended/visitor_stream.h"
#include "tools/gt_seqids.h"

static GtOptionParser* gt_seqids_option_parser_new(GT_UNUSED
                                                       void *tool_arguments)
{
  GtOptionParser *op;
  op = gt_option_parser_new("[GFF3_file]",
                            "Show sequence IDs from annotation file.");
  gt_option_parser_set_min_max_args(op, 1, 1);
  return op;
}

static int gt_seqids_runner(GT_UNUSED int argc, const char **argv,
                                  int parsed_args,
                                  GT_UNUSED void *tool_arguments, GtError *err)
{
  GtNodeStream *in_stream, *v_stream;
  GtCstrTable *cst;
  int had_err = 0;
  gt_error_check(err);

  cst = gt_cstr_table_new();
  in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                             argv + parsed_args);
  v_stream = gt_visitor_stream_new(in_stream, gt_collect_ids_visitor_new(cst));

  had_err = gt_node_stream_pull(v_stream, err);
  if (!had_err) {
    GtStrArray *seqids;
    unsigned long i;
    seqids = gt_cstr_table_get_all(cst);
    for (i = 0; i < gt_str_array_size(seqids); i++) {
      printf("%s\n", gt_str_array_get(seqids, i));
    }
    gt_str_array_delete(seqids);
  }

  gt_node_stream_delete(v_stream);
  gt_node_stream_delete(in_stream);
  gt_cstr_table_delete(cst);
  return had_err;
}

GtTool* gt_seqids(void)
{
  return gt_tool_new(NULL,
                  NULL,
                  gt_seqids_option_parser_new,
                  NULL,
                  gt_seqids_runner);
}
