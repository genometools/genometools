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

#include "core/io.h"
#include "core/ma.h"
#include "extended/bed_parser.h"

struct GtBEDParser {
  int placeholder;
};

GtBEDParser* gt_bed_parser_new(void)
{
  return gt_malloc(sizeof (GtBEDParser));
}

void gt_bed_parser_delete(GtBEDParser *bed_parser)
{
  if (!bed_parser) return;
  gt_free(bed_parser);
}

int gt_bed_parser_parse(GtBEDParser *bed_parser, GtQueue *genome_nodes,
                        const char *filename, GtError *err)
{
  GtIO *bed_file;
  gt_error_check(err);
  gt_assert(bed_parser && genome_nodes);
  bed_file = gt_io_new(filename, "r");
  gt_io_delete(bed_file);
  return 0;
}
