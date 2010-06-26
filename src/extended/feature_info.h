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

#ifndef FEATURE_INFO_H
#define FEATURE_INFO_H

#include "extended/genome_node.h"

typedef struct GtFeatureInfo GtFeatureInfo;

GtFeatureInfo* gt_feature_info_new(void);
void           gt_feature_info_delete(GtFeatureInfo*);
void           gt_feature_info_reset(GtFeatureInfo*);
GtFeatureNode* gt_feature_info_get(const GtFeatureInfo*, const char *id);
void           gt_feature_info_add(GtFeatureInfo*, const char *id,
                                   GtFeatureNode*);
GtFeatureNode* gt_feature_info_get_pseudo_parent(const GtFeatureInfo*,
                                                 const char *id);
void           gt_feature_info_add_pseudo_parent(GtFeatureInfo*, const char *id,
                                                 GtFeatureNode *pseudo_parent);
void           gt_feature_info_replace_pseudo_parent(GtFeatureInfo*,
                                                     GtFeatureNode *child,
                                                     GtFeatureNode
                                                     *new_pseudo_parent);
GtFeatureNode* gt_feature_info_find_root(const GtFeatureInfo*, const char *id);

#endif
