/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#ifndef FEATURE_INDEX_REP_H
#define FEATURE_INDEX_REP_H

#include <stdio.h>
#include "extended/feature_index.h"

typedef int         (*GtFeatureIndexAddRegionNodeFunc)(GtFeatureIndex*,
                                                       GtRegionNode*,
                                                       GtError*);
typedef int         (*GtFeatureIndexAddFeatureNodeFunc)(GtFeatureIndex*,
                                                        GtFeatureNode*,
                                                        GtError*);
typedef int         (*GtFeatureIndexRemoveNodeFunc)(GtFeatureIndex*,
                                                    GtFeatureNode*,
                                                    GtError*);
typedef GtArray*    (*GtFeatureIndexGetFeatsForSeqidFunc)(GtFeatureIndex*,
                                                          const char*,
                                                          GtError*);
typedef int         (*GtFeatureIndexGetFeatsForRangeFunc)(GtFeatureIndex*,
                                                          GtArray*,
                                                          const char*,
                                                          const GtRange*,
                                                          GtError*);
typedef char*       (*GtFeatureIndexGetFirstSeqidFunc)(const GtFeatureIndex*,
                                                       GtError*);
typedef int         (*GtFeatureIndexSaveFunc)(GtFeatureIndex*, GtError*);
typedef GtStrArray* (*GtFeatureIndexGetSeqidsFunc)(const GtFeatureIndex*,
                                                   GtError*);
typedef int         (*GtFeatureIndexGetRangeForSeqidFunc)(GtFeatureIndex*,
                                                          GtRange*,
                                                          const char*,
                                                          GtError*);
typedef int         (*GtFeatureIndexHasSeqidFunc)(const GtFeatureIndex*,
                                                  bool*,
                                                  const char*,
                                                  GtError*);
typedef void        (*GtFeatureIndexFreeFunc)(GtFeatureIndex*);

typedef struct GtFeatureIndexMembers GtFeatureIndexMembers;

struct GtFeatureIndex {
  const GtFeatureIndexClass *c_class;
  GtFeatureIndexMembers *pvt;
};

const GtFeatureIndexClass* gt_feature_index_class_new(size_t size,
                                         GtFeatureIndexAddRegionNodeFunc
                                                 add_region_node,
                                         GtFeatureIndexAddFeatureNodeFunc
                                                 add_feature_node,
                                         GtFeatureIndexRemoveNodeFunc
                                                 remove_node_func,
                                         GtFeatureIndexGetFeatsForSeqidFunc
                                                 get_features_for_seqid,
                                         GtFeatureIndexGetFeatsForRangeFunc
                                                 get_features_for_range,
                                         GtFeatureIndexGetFirstSeqidFunc
                                                 get_first_seqid,
                                         GtFeatureIndexSaveFunc
                                                 save_func,
                                         GtFeatureIndexGetSeqidsFunc
                                                 get_seqids,
                                         GtFeatureIndexGetRangeForSeqidFunc
                                                 get_range_for_seqid,
                                         GtFeatureIndexHasSeqidFunc
                                                 has_seqid,
                                         GtFeatureIndexFreeFunc
                                                 free);
GtFeatureIndex* gt_feature_index_create(const GtFeatureIndexClass*);
void*           gt_feature_index_cast(const GtFeatureIndexClass*,
                                      GtFeatureIndex*);

#endif
