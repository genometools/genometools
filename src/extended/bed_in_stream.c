/*
  Copyright (c) 2008-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/bed_in_stream.h"
#include "extended/bed_parser.h"
#include "extended/node_stream_api.h"

struct GtBEDInStream {
  const GtNodeStream parent_instance;
  GtBEDParser *bed_parser;
  GtQueue *genome_node_buffer;
  char *filename;
  bool file_processed;
};

#define bed_in_stream_cast(NS)\
        gt_node_stream_cast(gt_bed_in_stream_class(), NS)

static int bed_in_stream_process_file(GtBEDInStream *bed_in_stream,
                                      GtError *err)
{
  int had_err;
  gt_error_check(err);
  gt_assert(bed_in_stream);
  had_err = gt_bed_parser_parse(bed_in_stream->bed_parser,
                                bed_in_stream->genome_node_buffer,
                                bed_in_stream->filename, err);
  return had_err;
}

static int bed_in_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                              GT_UNUSED GtError *err)
{
  GtBEDInStream *is;
  int had_err = 0;
  gt_error_check(err);
  is = bed_in_stream_cast(ns);
  if (!is->file_processed) {
    had_err = bed_in_stream_process_file(is, err);
    is->file_processed = true;
  }
  if (!had_err && gt_queue_size(is->genome_node_buffer)) {
    /* we still have a node in the buffer -> serve it from there */
    *gn = gt_queue_get(is->genome_node_buffer);
    return 0;
  }
  if (!had_err) {
    /* the buffer is empty */
    gt_assert(!gt_queue_size(is->genome_node_buffer));
    *gn = NULL;
  }
  return had_err;
}

static void bed_in_stream_free(GtNodeStream *ns)
{
  GtBEDInStream *bed_in_stream = bed_in_stream_cast(ns);
  gt_free(bed_in_stream->filename);
  while (gt_queue_size(bed_in_stream->genome_node_buffer))
    gt_genome_node_delete(gt_queue_get(bed_in_stream->genome_node_buffer));
  gt_queue_delete(bed_in_stream->genome_node_buffer);
  gt_bed_parser_delete(bed_in_stream->bed_parser);
}

const GtNodeStreamClass* gt_bed_in_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtBEDInStream),
                                   bed_in_stream_free,
                                   bed_in_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_bed_in_stream_new(const char *filename)
{
  GtBEDInStream *bed_in_stream;
  GtNodeStream *ns = gt_node_stream_create(gt_bed_in_stream_class(), false);
  bed_in_stream = bed_in_stream_cast(ns);
  bed_in_stream->bed_parser = gt_bed_parser_new();
  bed_in_stream->genome_node_buffer = gt_queue_new();
  bed_in_stream->filename = filename ? gt_cstr_dup(filename) : NULL;
  return ns;
}

void gt_bed_in_stream_set_feature_type(GtBEDInStream *bed_in_stream,
                                       const char *type)
{
  gt_assert(bed_in_stream && type);
  gt_bed_parser_set_feature_type(bed_in_stream->bed_parser, type);
}

void gt_bed_in_stream_set_thick_feature_type(GtBEDInStream *bed_in_stream,
                                             const char *type)
{
  gt_assert(bed_in_stream && type);
  gt_bed_parser_set_thick_feature_type(bed_in_stream->bed_parser, type);
}

void gt_bed_in_stream_set_block_type(GtBEDInStream *bed_in_stream,
                                     const char *type)
{
  gt_assert(bed_in_stream && type);
  gt_bed_parser_set_block_type(bed_in_stream->bed_parser, type);
}
