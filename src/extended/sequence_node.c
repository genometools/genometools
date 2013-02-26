/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "core/unused_api.h"
#include "extended/sequence_node_api.h"
#include "extended/genome_node_rep.h"

struct GtSequenceNode {
  const GtGenomeNode parent_instance;
  GtStr *description,
        *sequence;
};

static void sequence_node_free(GtGenomeNode *gn)
{
  GtSequenceNode *sn = gt_sequence_node_cast(gn);
  gt_str_delete(sn->sequence);
  gt_str_delete(sn->description);
}

static GtStr* sequence_node_get_seqid(GtGenomeNode *gn)
{
  GtSequenceNode *sn = gt_sequence_node_cast(gn);
  return sn->description;
}

static GtRange sequence_node_get_range(GT_UNUSED GtGenomeNode *gn)
{
  GtRange range;
  range.start = 0;
  range.end = 0;
  return range;
}

static void sequence_node_change_seqid(GtGenomeNode *gn, GtStr *seqid)
{
  GtSequenceNode *sn = gt_sequence_node_cast(gn);
  gt_assert(sn && seqid);
  gt_str_delete(sn->description);
  sn->description = gt_str_ref(seqid);
}

static int sequence_node_accept(GtGenomeNode *gn, GtNodeVisitor *nv,
                                GtError *err)
{
  GtSequenceNode *sn;
  gt_error_check(err);
  sn = gt_sequence_node_cast(gn);
  return gt_node_visitor_visit_sequence_node(nv, sn, err);
}

const GtGenomeNodeClass* gt_sequence_node_class()
{
  static const GtGenomeNodeClass *gnc = NULL;
  gt_class_alloc_lock_enter();
  if (!gnc) {
    gnc = gt_genome_node_class_new(sizeof (GtSequenceNode),
                                   sequence_node_free,
                                   sequence_node_get_seqid,
                                   sequence_node_get_seqid,
                                   sequence_node_get_range,
                                   NULL,
                                   sequence_node_change_seqid,
                                   sequence_node_accept);
  }
  gt_class_alloc_lock_leave();
  return gnc;
}

GtGenomeNode* gt_sequence_node_new(const char *description, GtStr *sequence)
{
  GtGenomeNode *gn = gt_genome_node_create(gt_sequence_node_class());
  GtSequenceNode *sn = gt_sequence_node_cast(gn);
  gt_assert(description && sequence);
  sn->description = gt_str_new_cstr(description);
  sn->sequence = gt_str_ref(sequence);
  return gn;
}

const char* gt_sequence_node_get_description(const GtSequenceNode *sn)
{
  gt_assert(sn);
  return gt_str_get(sn->description);
}

const char* gt_sequence_node_get_sequence(const GtSequenceNode *sn)
{
  gt_assert(sn);
  return gt_str_get(sn->sequence);
}

unsigned long gt_sequence_node_get_sequence_length(const GtSequenceNode *sn)
{
  gt_assert(sn);
  return gt_str_length(sn->sequence);
}
