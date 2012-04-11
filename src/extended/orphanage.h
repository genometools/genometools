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

#ifndef ORPHANAGE_H
#define ORPHANAGE_H

typedef struct GtOrphanage GtOrphanage;

GtOrphanage*  gt_orphanage_new(void);
void          gt_orphanage_delete(GtOrphanage*);
void          gt_orphanage_reset(GtOrphanage*);
/* Takes ownership of <orphan*>. */
void          gt_orphanage_add(GtOrphanage*, GtGenomeNode *orphan,
                               const char *orphan_id,
                               GtStrArray *missing_parents);
void          gt_orphanage_reg_parent(GtOrphanage*, const char *id);
GtGenomeNode* gt_orphanage_get_orphan(GtOrphanage*);
bool          gt_orphanage_parent_is_missing(GtOrphanage*, const char *id);
bool          gt_orphanage_is_orphan(GtOrphanage*, const char *id);

#endif
