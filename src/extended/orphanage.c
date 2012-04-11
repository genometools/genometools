/*
  Copyright (c) 2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/cstr_api.h"
#include "core/cstr_table_api.h"
#include "core/hashmap_api.h"
#include "core/ma_api.h"
#include "core/queue_api.h"
#include "core/unused_api.h"
#include "extended/feature_node_api.h"
#include "extended/gff3_defines.h"
#include "extended/orphanage.h"

struct GtOrphanage {
  GtQueue *orphans;
  GtCstrTable *missing_parents,
              *orphan_ids;
};

GtOrphanage* gt_orphanage_new(void)
{
  GtOrphanage *o = gt_malloc(sizeof *o);
  o->orphans = gt_queue_new();
  o->missing_parents = gt_cstr_table_new();
  o->orphan_ids = gt_cstr_table_new();
  return o;
}

void gt_orphanage_delete(GtOrphanage *o)
{
  if (!o) return;
  gt_cstr_table_delete(o->orphan_ids);
  gt_cstr_table_delete(o->missing_parents);
  while (gt_queue_size(o->orphans))
    gt_genome_node_delete(gt_queue_get(o->orphans));
  gt_queue_delete(o->orphans);
  gt_free(o);
}

void gt_orphanage_reset(GtOrphanage *o)
{
  gt_assert(o);
  while (gt_queue_size(o->orphans))
    gt_genome_node_delete(gt_queue_get(o->orphans));
  gt_cstr_table_reset(o->missing_parents);
  gt_cstr_table_reset(o->orphan_ids);
}

void gt_orphanage_add(GtOrphanage *o, GtGenomeNode *orphan,
                      const char *orphan_id, GtStrArray *missing_parents)
{
  const char *missing_parent;
  unsigned long i;
  gt_assert(o && orphan);
  gt_assert(gt_feature_node_get_attribute((GtFeatureNode*) orphan,
                                          GT_GFF_PARENT));
  gt_queue_add(o->orphans, orphan);
  if (orphan_id && !gt_cstr_table_get(o->orphan_ids, orphan_id))
    gt_cstr_table_add(o->orphan_ids, orphan_id);
  if (missing_parents) {
    for (i = 0; i < gt_str_array_size(missing_parents); i++) {
      missing_parent = gt_str_array_get(missing_parents, i);
      if (!gt_cstr_table_get(o->missing_parents, missing_parent))
        gt_cstr_table_add(o->missing_parents, missing_parent);
    }
  }
}

void gt_orphanage_reg_parent(GtOrphanage *o, const char *parent_id)
{
  gt_assert(o && parent_id);
  gt_cstr_table_remove(o->missing_parents, parent_id);
}

GtGenomeNode* gt_orphanage_get_orphan(GtOrphanage *o)
{
  gt_assert(o);
  if (gt_queue_size(o->orphans))
    return gt_queue_get(o->orphans);
  return NULL;
}

bool gt_orphanage_parent_is_missing(GtOrphanage *o, const char *parent_id)
{
  gt_assert(o && parent_id);
  if (gt_cstr_table_get(o->missing_parents, parent_id))
    return true;
  return false;
}

bool gt_orphanage_is_orphan(GtOrphanage *o, const char *id)
{
 gt_assert(o && id);
  if (gt_cstr_table_get(o->orphan_ids, id))
    return true;
  return false;
}
