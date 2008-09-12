/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
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
#include "core/ma.h"
#include "core/unused_api.h"
#include "annotationsketch/feature_index_rep.h"

GtFeatureIndex* gt_feature_index_create(const GtFeatureIndexClass *fic)
{
  GtFeatureIndex *fi;
  assert(fic && fic->size);
  fi = gt_calloc(1, fic->size);
  fi->c_class = fic;
  return fi;
}

GtFeatureIndex* gt_feature_index_ref(GtFeatureIndex *fi)
{
  assert(fi);
  fi->reference_count++;
  return fi;
}

void gt_feature_index_delete(GtFeatureIndex *fi)
{
  if (!fi) return;
  if (fi->reference_count) {
    fi->reference_count--;
    return;
  }
  assert(fi->c_class);
  if (fi->c_class->free)
    fi->c_class->free(fi);
  gt_free(fi);
}

void gt_feature_index_add_region_node(GtFeatureIndex *fi, GtRegionNode *rn)
{
  assert(fi && fi->c_class && rn);
  fi->c_class->add_region_node(fi, rn);
}

void gt_feature_index_add_genome_feature(GtFeatureIndex *fi,
                                         GtFeatureNode *fn)
{
  assert(fi && fi->c_class && fn);
  fi->c_class->add_genome_feature(fi, fn);
}

int gt_feature_index_add_gff3file(GtFeatureIndex *fi,
                                  const char *gff3file,
                                  GtError *err)
{
  assert(fi && fi->c_class && gff3file);
  return fi->c_class->add_gff3file(fi, gff3file, err);
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
                                            GtRange rng , GtError *err)
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

void* gt_feature_index_cast(const GtFeatureIndexClass *fic,
                            GtFeatureIndex *fi)
{
  assert(fic && fi && fi->c_class == fic);
  return fi;
}
