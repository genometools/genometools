/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/array.h"
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/unused_api.h"
#include "extended/chseqids_stream.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "extended/mapping.h"
#include "extended/region_node.h"

struct GtChseqidsStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtMapping *chseqids_mapping;
  GtArray *node_buffer;
  unsigned long buffer_index;
  bool sequence_regions_processed;
};

#define chseqids_stream_cast(GS)\
        gt_node_stream_cast(gt_chseqids_stream_class(), GS)

static int change_sequence_id(GtFeatureNode *fn, void *data,
                              GT_UNUSED GtError *err)
{
  GtStr *changed_seqid = data;
  gt_error_check(err);
  gt_assert(changed_seqid);
  gt_genome_node_change_seqid((GtGenomeNode*) fn, changed_seqid);
  return 0;
}

static int chseqids_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                GtError *err)
{
  GtChseqidsStream *cs;
  GtGenomeNode *node, **gn_a, **gn_b;
  GtFeatureNode *feature_node;
  GtStr *changed_seqid;
  unsigned long i;
  int GT_UNUSED rval, had_err = 0;
  gt_error_check(err);
  cs = chseqids_stream_cast(ns);

  if (!cs->sequence_regions_processed) {
    while (!had_err) {
      if (!(had_err = gt_node_stream_next(cs->in_stream, &node, err))) {
        if (node)
          gt_array_add(cs->node_buffer, node);
        else
          break;
        if (!gt_region_node_try_cast(node))
          break; /* no more sequence regions */
      }
    }
    /* now the buffer contains only sequence regions (except the last entry)
       -> change sequence ids */
    for (i = 0; !had_err && i < gt_array_size(cs->node_buffer); i++) {
      node = *(GtGenomeNode**) gt_array_get(cs->node_buffer, i);
      if (gt_genome_node_get_seqid(node)) {
        if  ((changed_seqid = gt_mapping_map_string(cs->chseqids_mapping,
                                     gt_str_get(gt_genome_node_get_seqid(node)),
                                                 err))) {
          if ((feature_node = gt_feature_node_try_cast(node))) {
            rval = gt_feature_node_traverse_children(feature_node,
                                                     changed_seqid,
                                                     change_sequence_id, true,
                                                     err);
            gt_assert(!rval); /* change_sequence_id() is sane */
          }
          else
            gt_genome_node_change_seqid(node, changed_seqid);
          gt_str_delete(changed_seqid);
        }
        else
          had_err = -1;
       }
    }
    /* sort them */
    if (!had_err)
      gt_genome_nodes_sort(cs->node_buffer);
    /* consolidate them */
    for (i = 1; !had_err && i + 1 < gt_array_size(cs->node_buffer); i++) {
      gn_a = gt_array_get(cs->node_buffer, i-1);
      gn_b = gt_array_get(cs->node_buffer, i);
      if (gt_genome_nodes_are_equal_region_nodes(*gn_a, *gn_b)) {
        gt_region_node_consolidate(gt_region_node_cast(*gn_b),
                                   gt_region_node_cast(*gn_a));
        gt_genome_node_delete(*gn_a);
        *gn_a = NULL;
      }
    }
    cs->sequence_regions_processed = true;
  }

  /* return non-null nodes from buffer */
  while (!had_err && cs->buffer_index < gt_array_size(cs->node_buffer)) {
    node = *(GtGenomeNode**) gt_array_get(cs->node_buffer, cs->buffer_index);
    cs->buffer_index++;
    if (node) {
      *gn = node;
      return had_err;
    }
  }

  if (!had_err)
    had_err = gt_node_stream_next(cs->in_stream, gn, err);
  if (!had_err && *gn) {
    if (gt_genome_node_get_seqid(*gn)) {
      changed_seqid = gt_mapping_map_string(cs->chseqids_mapping,
                                      gt_str_get(gt_genome_node_get_seqid(*gn)),
                                         err);
      gt_assert(changed_seqid); /* is always defined, because an undefined
                                   mapping would be catched earlier */
      if ((feature_node = gt_feature_node_try_cast(*gn))) {
        rval = gt_feature_node_traverse_children(feature_node, changed_seqid,
                                                 change_sequence_id, true, err);
        gt_assert(!rval); /* change_sequence_id() is sane */
      }
      else
        gt_genome_node_change_seqid(*gn, changed_seqid);
      gt_str_delete(changed_seqid);
    }
  }

  return had_err;
}

static void chseqids_stream_free(GtNodeStream *ns)
{
  GtChseqidsStream *cs;
  unsigned long i;
  cs = chseqids_stream_cast(ns);
  gt_mapping_delete(cs->chseqids_mapping);
  for (i = cs->buffer_index; i < gt_array_size(cs->node_buffer); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(cs->node_buffer, i));
  gt_array_delete(cs->node_buffer);
  gt_node_stream_delete(cs->in_stream);
}

const GtNodeStreamClass* gt_chseqids_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtChseqidsStream),
                                   chseqids_stream_free,
                                   chseqids_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_chseqids_stream_new(GtNodeStream *in_stream,
                                     GtStr *chseqids_file, GtError *err)
{
  GtNodeStream *ns;
  GtChseqidsStream *cs;
  gt_error_check(err);
  gt_assert(in_stream && chseqids_file);
  gt_assert(gt_node_stream_is_sorted(in_stream));
  ns = gt_node_stream_create(gt_chseqids_stream_class(), false);
  cs = chseqids_stream_cast(ns);
  cs->in_stream = gt_node_stream_ref(in_stream);
  cs->chseqids_mapping = gt_mapping_new(chseqids_file, "chseqids",
                                        GT_MAPPINGTYPE_STRING, err);
  if (!cs->chseqids_mapping) {
    gt_node_stream_delete(ns);
    return NULL;
  }
  cs->node_buffer = gt_array_new(sizeof (GtGenomeNode*));
  return ns;
}
