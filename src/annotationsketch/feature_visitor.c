/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <assert.h>
#include "core/unused_api.h"
#include "extended/genome_visitor_rep.h"
#include "extended/region_node.h"
#include "annotationsketch/feature_index.h"
#include "annotationsketch/feature_visitor.h"

struct FeatureVisitor {
  const GenomeVisitor parent_instance;
        GtFeatureIndex *feature_index;
};

#define feature_visitor_cast(GV)\
        genome_visitor_cast(feature_visitor_class(), GV)

static void feature_visitor_free(GenomeVisitor *gv)
{
  FeatureVisitor *feature_visitor = feature_visitor_cast(gv);
  assert(feature_visitor);
  gt_feature_index_delete(feature_visitor->feature_index);
}

static int feature_visitor_genome_feature(GenomeVisitor *gv,
                                          GtGenomeFeature *gf,
                                          GT_UNUSED GtError *err)
{
  FeatureVisitor *v = feature_visitor_cast(gv);
  gt_error_check(err);
  gt_feature_index_add_genome_feature(v->feature_index, gf);
  return 0;
}

static int feature_visitor_region_node(GenomeVisitor *gv, GtRegionNode *rn,
                                       GT_UNUSED GtError *err)
{
  FeatureVisitor *v = feature_visitor_cast(gv);
  gt_error_check(err);
  gt_feature_index_add_region_node(v->feature_index, rn);
  return 0;
}

const GenomeVisitorClass* feature_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (FeatureVisitor),
                                          feature_visitor_free,
                                          NULL,
                                          feature_visitor_genome_feature,
                                          feature_visitor_region_node,
                                          NULL };
  return &gvc;
}

GenomeVisitor* feature_visitor_new(GtFeatureIndex *fi)
{
  GenomeVisitor *gv;
  FeatureVisitor *feature_visitor;
  assert(fi != NULL);
  gv = genome_visitor_create(feature_visitor_class());
  feature_visitor = feature_visitor_cast(gv);
  feature_visitor->feature_index = gt_feature_index_ref(fi);
  assert(feature_visitor != NULL);
  return gv;
}
