/*
  Copyright (c) 2006-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/fa.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/gtf_in_stream.h"
#include "extended/gtf_parser.h"
#include "extended/node_stream_api.h"
#include "extended/type_checker_builtin_api.h"

struct GtGTFInStream {
  const GtNodeStream parent_instance;
  GtQueue *genome_node_buffer;
  GtTypeChecker *type_checker;
  char *filename;
  bool file_processed,
       tidy;
};

#define gtf_in_stream_cast(NS)\
        gt_node_stream_cast(gt_gtf_in_stream_class(), NS)

static int gtf_in_stream_process_file(GtGTFInStream *gtf_in_stream,
                                      GtError *err)
{
  GtGTFParser *gtf_parser;
  GtStr *filenamestr;
  GtFile *fpin;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(gtf_in_stream);

  gtf_parser = gt_gtf_parser_new(gtf_in_stream->type_checker);

  /* open input file */
  if (gtf_in_stream->filename) {
    if (!(fpin = gt_file_new(gtf_in_stream->filename, "r", err)))
      had_err = -1;
  }
  else
    fpin = NULL;

  /* parse input file */
  if (!had_err) {
    filenamestr = gt_str_new_cstr(gtf_in_stream->filename
                                  ? gtf_in_stream->filename : "stdin");
    had_err = gt_gtf_parser_parse(gtf_parser, gtf_in_stream->genome_node_buffer,
                                  filenamestr, fpin, gtf_in_stream->tidy, err);
    gt_str_delete(filenamestr);
  }

  /* close input file, if necessary */
  gt_file_delete(fpin);

  /* free */
  gt_gtf_parser_delete(gtf_parser);

  return had_err;
}

static int gtf_in_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                              GT_UNUSED GtError *err)
{
  GtGTFInStream *is;
  int had_err = 0;
  gt_error_check(err);
  is = gtf_in_stream_cast(ns);
  if (!is->file_processed) {
    had_err = gtf_in_stream_process_file(is, err);
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

static void gtf_in_stream_free(GtNodeStream *ns)
{
  GtGTFInStream *gtf_in_stream = gtf_in_stream_cast(ns);
  gt_free(gtf_in_stream->filename);
  gt_type_checker_delete(gtf_in_stream->type_checker);
  while (gt_queue_size(gtf_in_stream->genome_node_buffer))
    gt_genome_node_delete(gt_queue_get(gtf_in_stream->genome_node_buffer));
  gt_queue_delete(gtf_in_stream->genome_node_buffer);
}

const GtNodeStreamClass* gt_gtf_in_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtGTFInStream),
                                   gtf_in_stream_free,
                                   gtf_in_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_gtf_in_stream_new(const char *filename)
{
  GtGTFInStream *gtf_in_stream;
  GtNodeStream *ns = gt_node_stream_create(gt_gtf_in_stream_class(), false);
  gtf_in_stream = gtf_in_stream_cast(ns);
  gtf_in_stream->genome_node_buffer = gt_queue_new();
  gtf_in_stream->type_checker = gt_type_checker_builtin_new();
  gtf_in_stream->filename = filename ? gt_cstr_dup(filename) : NULL;
  return ns;
}

void gt_gtf_in_stream_enable_tidy_mode(GtNodeStream *ns)
{
  GtGTFInStream *is = gtf_in_stream_cast(ns);
  is->tidy = true;
}
