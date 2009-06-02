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

#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/array.h"
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
  gt_assert(fic && fic->size);
  fi = gt_calloc(1, fic->size);
  fi->c_class = fic;
  fi->pvt = gt_calloc(1, sizeof (GtFeatureIndexMembers));
  return fi;
}

GtFeatureIndex* gt_feature_index_ref(GtFeatureIndex *fi)
{
  gt_assert(fi);
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
  gt_assert(fi->c_class);
  if (fi->c_class->free)
    fi->c_class->free(fi);
  gt_free(fi->pvt);
  gt_free(fi);
}

void gt_feature_index_add_region_node(GtFeatureIndex *feature_index,
                                      GtRegionNode *region_node)
{
  gt_assert(feature_index && feature_index->c_class && region_node);
  feature_index->c_class->add_region_node(feature_index, region_node);
}

void gt_feature_index_add_feature_node(GtFeatureIndex *feature_index,
                                       GtFeatureNode *feature_node)
{
  gt_assert(feature_index && feature_index->c_class && feature_node);
  feature_index->c_class->add_feature_node(feature_index, feature_node);
}

int gt_feature_index_add_gff3file(GtFeatureIndex *feature_index,
                                  const char *gff3file, GtError *err)
{
  GtNodeStream *gff3_in_stream;
  GtGenomeNode *gn;
  GtArray *tmp;
  int had_err = 0;
  unsigned long i;
  gt_error_check(err);
  gt_assert(feature_index && gff3file);
  tmp = gt_array_new(sizeof (GtGenomeNode*));
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(1, &gff3file);
  while (!(had_err = gt_node_stream_next(gff3_in_stream, &gn, err)) && gn)
    gt_array_add(tmp, gn);
  if (!had_err) {
    GtNodeVisitor  *feature_visitor = gt_feature_visitor_new(feature_index);
    for (i=0;i<gt_array_size(tmp);i++) {
      gn = *(GtGenomeNode**) gt_array_get(tmp, i);
      had_err = gt_genome_node_accept(gn, feature_visitor, NULL);
      gt_assert(!had_err); /* cannot happen */
    }
    gt_node_visitor_delete(feature_visitor);
  }
  gt_node_stream_delete(gff3_in_stream);
  for (i=0;i<gt_array_size(tmp);i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(tmp, i));
  gt_array_delete(tmp);
  return had_err;
}

GtArray* gt_feature_index_get_features_for_seqid(GtFeatureIndex *fi,
                                                 const char *seqid)
{
  gt_assert(fi && fi->c_class && seqid);
  return fi->c_class->get_features_for_seqid(fi, seqid);
}

int gt_feature_index_get_features_for_range(GtFeatureIndex *feature_index,
                                            GtArray *results,
                                            const char *seqid,
                                            const GtRange *range, GtError *err)
{
  gt_assert(feature_index && feature_index->c_class && results && seqid &&
            range);
  gt_assert(gt_range_length(range) > 0);
  return feature_index->c_class->get_features_for_range(feature_index, results,
                                                        seqid, range, err);
}

const char* gt_feature_index_get_first_seqid(const GtFeatureIndex
                                              *feature_index)
{
  gt_assert(feature_index && feature_index->c_class);
  return feature_index->c_class->get_first_seqid(feature_index);
}

GtStrArray* gt_feature_index_get_seqids(const GtFeatureIndex *feature_index)
{
  gt_assert(feature_index && feature_index->c_class);
  return feature_index->c_class->get_seqids(feature_index);
}

void gt_feature_index_get_range_for_seqid(GtFeatureIndex *feature_index,
                                          GtRange *range, const char *seqid)
{
  gt_assert(feature_index && feature_index->c_class && range && seqid);
  feature_index->c_class->get_range_for_seqid(feature_index, range, seqid);
}

bool gt_feature_index_has_seqid(const GtFeatureIndex *feature_index,
                                const char *seqid)
{
  gt_assert(feature_index && feature_index->c_class && seqid);
  return feature_index->c_class->has_seqid(feature_index, seqid);
}

void* gt_feature_index_cast(GT_UNUSED const GtFeatureIndexClass *fic,
                            GtFeatureIndex *fi)
{
  gt_assert(fic && fi && fi->c_class == fic);
  return fi;
}
