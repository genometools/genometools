/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef FEATURE_NODE_OBSERVER_H
#define FEATURE_NODE_OBSERVER_H

#include "core/phase_api.h"
#include "core/str_api.h"
#include "core/strand_api.h"

typedef struct GtFeatureNodeObserver GtFeatureNodeObserver;

#include "extended/feature_node.h"

typedef void (*GtFeatureNodeScoreChangedCallback)(GtFeatureNode *observed,
                                                  double score,
                                                  void *data);
typedef void (*GtFeatureNodePhaseChangedCallback)(GtFeatureNode *observed,
                                                  GtPhase phase,
                                                  void *data);
typedef void (*GtFeatureNodeStrandChangedCallback)(GtFeatureNode *observed,
                                                   GtStrand strand,
                                                   void *data);
typedef void (*GtFeatureNodeTypeChangedCallback)(GtFeatureNode *observed,
                                                 const char *type,
                                                 void *data);
typedef void (*GtFeatureNodeSourceChangedCallback)(GtFeatureNode *observed,
                                                   GtStr *source,
                                                   void *data);
typedef void (*GtFeatureNodeAttributeChangedCallback)(GtFeatureNode *observed,
                                                      bool added,
                                                      const char *key,
                                                      const char *value,
                                                      void *data);
typedef void (*GtFeatureNodeAttributeDeletedCallback)(GtFeatureNode *observed,
                                                      const char *key,
                                                      void *data);
typedef void (*GtFeatureNodeRangeChangedCallback)(GtFeatureNode *observed,
                                                  GtRange *range,
                                                  void *data);
typedef void (*GtFeatureNodeMultiChangedCallback)(GtFeatureNode *observed,
                                                  bool is_multi,
                                                  GtFeatureNode *repr,
                                                  void *data);
typedef void (*GtFeatureNodeMarkChangedCallback)(GtFeatureNode *observed,
                                                 bool marked,
                                                 void *data);
typedef int  (*GtFeatureNodeChildAddedCallback)(GtFeatureNode *observed,
                                                GtFeatureNode *child,
                                                void *data);
typedef void (*GtFeatureNodeDeletedCallback)(GtFeatureNode *observed,
                                             void *data);

struct GtFeatureNodeObserver {
  GtFeatureNodeScoreChangedCallback score_changed;
  GtFeatureNodePhaseChangedCallback phase_changed;
  GtFeatureNodeStrandChangedCallback strand_changed;
  GtFeatureNodeTypeChangedCallback type_changed;
  GtFeatureNodeSourceChangedCallback source_changed;
  GtFeatureNodeAttributeChangedCallback attribute_changed;
  GtFeatureNodeAttributeDeletedCallback attribute_deleted;
  GtFeatureNodeRangeChangedCallback range_changed;
  GtFeatureNodeMultiChangedCallback multi_changed;
  GtFeatureNodeMarkChangedCallback mark_changed;
  GtFeatureNodeChildAddedCallback child_added;
  GtFeatureNodeDeletedCallback deleted;
  void *data;
  unsigned long reference_count;
};

GtFeatureNodeObserver* gt_feature_node_observer_new(void);
GtFeatureNodeObserver* gt_feature_node_observer_ref(GtFeatureNodeObserver *o);
void                   gt_feature_node_observer_delete(GtFeatureNodeObserver
                                                                            *o);

#endif
