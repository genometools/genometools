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
#include "core/cstr_table_api.h"
#include "core/error_api.h"
#include "core/ma.h"
#include "core/queue_api.h"
#include "core/str_array_api.h"
#include "core/unused_api.h"
#include "extended/node_stream_api.h"
#include "extended/collect_ids_visitor.h"
#include "extended/node_stream_api.h"
#include "extended/region_mapping_api.h"
#include "extended/sequence_node_add_stream.h"

struct GtSequenceNodeAddStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *collect_vis;
  GtCstrTable *seqid_table;
  GtStrArray *seqids;
  GtRegionMapping *rm;
  GtUword cur_seqid;
};

const GtNodeStreamClass* gt_sequence_node_add_stream_class(void);

#define gt_sequence_node_add_stream_cast(GS)\
        gt_node_stream_cast(gt_sequence_node_add_stream_class(), GS);

static int sequence_node_add_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                         GtError *err)
{
  GtSequenceNodeAddStream *s;
  int had_err;
  gt_error_check(err);
  s = gt_sequence_node_add_stream_cast(ns);

  /* stream nodes as long as we have some, record seen seqids */
  if (!(had_err = gt_node_stream_next(s->in_stream, gn, err)) && *gn) {
    had_err = gt_genome_node_accept(*gn, s->collect_vis, err);
  }

  /* if there are no more  */
  if (!had_err && !*gn) {
    if (!s->seqids) {
      s->seqids = gt_cstr_table_get_all(s->seqid_table);
    }
    gt_assert(s->seqids);
    if (s->cur_seqid >= gt_str_array_size(s->seqids)) {
      *gn = NULL;
      return 0;
    } else {
      GtGenomeNode *new_sn;
      GtUword len;
      char *seq = NULL;
      GtStr *seqid = gt_str_new(),
            *seqstr = gt_str_new();
      gt_str_append_cstr(seqid, gt_str_array_get(s->seqids, s->cur_seqid));
      had_err = gt_region_mapping_get_sequence_length(s->rm, &len, seqid, err);
      if (!had_err) {

        had_err = gt_region_mapping_get_sequence(s->rm, &seq, seqid, 1, len,
                                                 err);
      }
      if (!had_err) {
        gt_str_append_cstr(seqstr, seq);
        new_sn = gt_sequence_node_new(gt_str_get(seqid), seqstr);
        *gn = new_sn;
      }
      s->cur_seqid++;
      gt_free(seq);
      gt_str_delete(seqid);
      gt_str_delete(seqstr);
    }
  }

  return had_err;
}

static void sequence_node_add_stream_free(GtNodeStream *ns)
{
  GtSequenceNodeAddStream *fs = gt_sequence_node_add_stream_cast(ns);
  gt_node_visitor_delete(fs->collect_vis);
  gt_cstr_table_delete(fs->seqid_table);
  gt_str_array_delete(fs->seqids);
  gt_node_stream_delete(fs->in_stream);
}

const GtNodeStreamClass* gt_sequence_node_add_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtSequenceNodeAddStream),
                                   sequence_node_add_stream_free,
                                   sequence_node_add_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_sequence_node_add_stream_new(GtNodeStream *in_stream,
                                              GtRegionMapping *rm,
                                              GT_UNUSED GtError *err)
{
  GtNodeStream *ns = gt_node_stream_create(gt_sequence_node_add_stream_class(),
                                           gt_node_stream_is_sorted(in_stream));
  GtSequenceNodeAddStream *s = gt_sequence_node_add_stream_cast(ns);
  gt_assert(in_stream);
  s->rm = rm;
  s->in_stream = gt_node_stream_ref(in_stream);
  s->seqid_table = gt_cstr_table_new();
  s->seqids = NULL;
  s->cur_seqid = 0;
  s->collect_vis = gt_collect_ids_visitor_new(s->seqid_table);
  return ns;
}
