/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd

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
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/node_stream_api.h"
#include "extended/sequence_node_out_stream.h"
#include "extended/sequence_node_out_visitor.h"
#include "core/error_api.h"

struct GtSequenceNodeOutStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *sequence_node_out_visitor;
};

const GtNodeStreamClass* gt_sequence_node_out_stream_class(void);

#define gt_sequence_node_out_stream_cast(GS)\
        gt_node_stream_cast(gt_sequence_node_out_stream_class(), GS);

static int sequence_node_out_stream_next(GtNodeStream *ns, GtGenomeNode **gn, GtError *err)
{
  GtSequenceNodeOutStream *fs;
  int had_err;
  gt_error_check(err);
  fs = gt_sequence_node_out_stream_cast(ns);

  /* we still have nodes in the buffer */
  if (gt_sequence_node_out_visitor_node_buffer_size(fs->sequence_node_out_visitor)) {
    /* return one of them */
    *gn = gt_sequence_node_out_visitor_get_node(fs->sequence_node_out_visitor);
    return 0;
  }

  /* no nodes in the buffer -> get new nodes */
  while (!(had_err = gt_node_stream_next(fs->in_stream, gn, err)) && *gn) {
    gt_assert(*gn && !had_err);
    had_err = gt_genome_node_accept(*gn, fs->sequence_node_out_visitor, err);
    if (had_err) {
      /* we own the node -> delete it */
      gt_genome_node_delete(*gn);
      *gn = NULL;
      break;
    }
    if (gt_sequence_node_out_visitor_node_buffer_size(fs->sequence_node_out_visitor)) {
      *gn = gt_sequence_node_out_visitor_get_node(fs->sequence_node_out_visitor);
      return 0;
    }
  }

  /* either we have an error or no new node */
  gt_assert(had_err || !*gn);
  return had_err;
}

static void sequence_node_out_stream_free(GtNodeStream *ns)
{
  GtSequenceNodeOutStream *fs = gt_sequence_node_out_stream_cast(ns);
  gt_node_visitor_delete(fs->sequence_node_out_visitor);
  gt_node_stream_delete(fs->in_stream);
}

const GtNodeStreamClass* gt_sequence_node_out_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtSequenceNodeOutStream),
                                   sequence_node_out_stream_free,
                                   sequence_node_out_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_sequence_node_out_stream_new(GtNodeStream *in_stream,
                                              GtFile *outfile,
                                              GT_UNUSED GtError *err)
{
  GtNodeStream *ns = gt_node_stream_create(gt_sequence_node_out_stream_class(),
                                           gt_node_stream_is_sorted(in_stream));
  GtSequenceNodeOutStream *s = gt_sequence_node_out_stream_cast(ns);
  gt_assert(in_stream);
  s->in_stream = gt_node_stream_ref(in_stream);
  s->sequence_node_out_visitor =
    gt_sequence_node_out_visitor_new(outfile);
  if  (!s->sequence_node_out_visitor) {
    gt_node_stream_delete(ns);
    return NULL;
  }
  return ns;
}
