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

#include <assert.h>
#include "core/class_alloc.h"
#include "core/queue.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "annotationsketch/feature_index_rep.h"
#include "annotationsketch/feature_visitor.h"
#include "extended/gff3_in_stream.h"

struct GtFeatureIndexClass {
  size_t size;
  GtFeatureIndexAddRegionNodeFunc add_region_node;
  GtFeatureIndexAddFeatureNodeFunc add_feature_node;
  GtFeatureIndexGetFeatsForSeqidFunc get_features_for_seqid;
  GtFeatureIndexGetFeatsForRangeFunc get_features_for_range;
  GtFeatureIndexGetFirstSeqidFunc get_first_seqid;
  GtFeatureIndexGetSeqidsFunc get_seqids;
  GtFeatureIndexGetRangeForSeqidFunc get_range_for_seqid;
  GtFeatureIndexHasSeqidFunc has_seqid;
  GtFeatureIndexFreeFunc free;
};

struct GtFeatureIndexMembers {
  unsigned int reference_count;
};

const GtFeatureIndexClass* gt_feature_index_class_new(size_t size,
                                         GtFeatureIndexAddRegionNodeFunc
                                                 add_region_node,
                                         GtFeatureIndexAddFeatureNodeFunc
                                                 add_feature_node,
                                         GtFeatureIndexGetFeatsForSeqidFunc
                                                 get_features_for_seqid,
                                         GtFeatureIndexGetFeatsForRangeFunc
                                                 get_features_for_range,
                                         GtFeatureIndexGetFirstSeqidFunc
                                                 get_first_seqid,
                                         GtFeatureIndexGetSeqidsFunc
                                                 get_seqids,
                                         GtFeatureIndexGetRangeForSeqidFunc
                                                 get_range_for_seqid,
                                         GtFeatureIndexHasSeqidFunc
                                                 has_seqid,
                                         GtFeatureIndexFreeFunc
                                                 free)
{
  GtFeatureIndexClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->add_region_node = add_region_node;
  c_class->add_feature_node = add_feature_node;
  c_class->get_features_for_seqid = get_features_for_seqid;
  c_class->get_features_for_range = get_features_for_range;
  c_class->get_first_seqid = get_first_seqid;
  c_class->get_seqids = get_seqids;
  c_class->get_range_for_seqid = get_range_for_seqid;
  c_class->has_seqid = has_seqid;
  c_class->free = free;
  return c_class;
}

GtFeatureIndex* gt_feature_index_create(const GtFeatureIndexClass *fic)
{
  GtFeatureIndex *fi;
  assert(fic && fic->size);
  fi = gt_calloc(1, fic->size);
  fi->c_class = fic;
  fi->pvt = gt_calloc(1, sizeof (GtFeatureIndexMembers));
  return fi;
}

GtFeatureIndex* gt_feature_index_ref(GtFeatureIndex *fi)
{
  assert(fi);
  fi->pvt->reference_count++;
  return fi;
}

void gt_feature_index_delete(GtFeatureIndex *fi)
{
  if (!fi) return;
  if (fi->pvt->reference_count) {
    fi->pvt->reference_count--;
    return;
  }
  assert(fi->c_class);
  if (fi->c_class->free)
    fi->c_class->free(fi);
  gt_free(fi->pvt);
  gt_free(fi);
}

void gt_feature_index_add_region_node(GtFeatureIndex *fi, GtRegionNode *rn)
{
  assert(fi && fi->c_class && rn);
  fi->c_class->add_region_node(fi, rn);
}

void gt_feature_index_add_feature_node(GtFeatureIndex *fi, GtFeatureNode *fn)
{
  assert(fi && fi->c_class && fn);
  fi->c_class->add_feature_node(fi, fn);
}

int gt_feature_index_add_gff3file(GtFeatureIndex *fi,
                                  const char *gff3file, GtError *err)
{
  GtNodeStream *gff3_in_stream;
  GtGenomeNode *gn;
  GtQueue *queue;
  int had_err = 0;
  gt_error_check(err);
  assert(fi && gff3file);
  queue = gt_queue_new();
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(1, &gff3file, false, false);
  while (!(had_err = gt_node_stream_next(gff3_in_stream, &gn, err)) && gn)
    gt_queue_add(queue, gn);
  if (!had_err) {
    GtNodeVisitor  *feature_visitor = gt_feature_visitor_new(fi);
    while (gt_queue_size(queue)) {
      gn = gt_queue_get(queue);
      had_err = gt_genome_node_accept(gn, feature_visitor, NULL);
      assert(!had_err); /* cannot happen */
    }
    gt_node_visitor_delete(feature_visitor);
  }
  gt_node_stream_delete(gff3_in_stream);
  while (gt_queue_size(queue))
    gt_genome_node_rec_delete(gt_queue_get(queue));
  gt_queue_delete(queue);
  return had_err;
}

GtArray* gt_feature_index_get_features_for_seqid(GtFeatureIndex *fi,
                                                 const char *seqid)
{
  assert(fi && fi->c_class && seqid);
  return fi->c_class->get_features_for_seqid(fi, seqid);
}

int gt_feature_index_get_features_for_range(GtFeatureIndex *fi,
                                            GtArray *results,
                                            const char *seqid,
                                            const GtRange *rng, GtError *err)
{
  assert(fi && fi->c_class && results && seqid);
  return fi->c_class->get_features_for_range(fi, results, seqid, rng, err);
}

const char* gt_feature_index_get_first_seqid(const GtFeatureIndex *fi)
{
  assert(fi && fi->c_class);
  return fi->c_class->get_first_seqid(fi);
}

GtStrArray* gt_feature_index_get_seqids(const GtFeatureIndex *fi)
{
  assert(fi && fi->c_class);
  return fi->c_class->get_seqids(fi);
}

void gt_feature_index_get_range_for_seqid(GtFeatureIndex *fi, GtRange *rng,
                                          const char *seqid)
{
  assert(fi && fi->c_class && rng && seqid);
  fi->c_class->get_range_for_seqid(fi, rng, seqid);
}

bool gt_feature_index_has_seqid(const GtFeatureIndex *fi, const char *seqid)
{
  assert(fi && fi->c_class && seqid);
  return fi->c_class->has_seqid(fi, seqid);
}

void* gt_feature_index_cast(GT_UNUSED const GtFeatureIndexClass *fic,
                            GtFeatureIndex *fi)
{
  assert(fic && fi && fi->c_class == fic);
  return fi;
}
