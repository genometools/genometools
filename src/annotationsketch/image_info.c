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

#include "annotationsketch/image_info.h"
#include "core/array.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/range.h"
#include "core/strand_api.h"
#include "core/unused_api.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"

struct GtImageInfo {
  GtArray* recmaps;
  unsigned int height;
  unsigned long reference_count;
};

GtImageInfo* gt_image_info_new()
{
  GtImageInfo *ii;
  ii = gt_calloc(1, sizeof (GtImageInfo));
  ii->recmaps = gt_array_new(sizeof (GtRecMap*));
  gt_assert(ii->recmaps);
  return ii;
}

GtImageInfo* gt_image_info_ref(GtImageInfo *ii)
{
  gt_assert(ii);
  ii->reference_count++;
  return ii;
}

void gt_image_info_delete(GtImageInfo *image_info)
{
  unsigned long i;
  if (!image_info) return;
  if (image_info->reference_count) {
    image_info->reference_count--;
    return;
  }
  for (i=0;i<gt_image_info_num_of_rec_maps(image_info);i++)
  {
    GtRecMap *rm = *(GtRecMap**) gt_array_get(image_info->recmaps, i);
    gt_rec_map_delete(rm);
  }
  gt_array_delete(image_info->recmaps);
  gt_free(image_info);
}

void gt_image_info_add_rec_map(GtImageInfo *ii, GtRecMap *rec_map)
{
  gt_assert(ii && rec_map);
  gt_array_add(ii->recmaps, rec_map);
}

void gt_image_info_set_height(GtImageInfo *ii, unsigned int height)
{
  gt_assert(ii);
  ii->height = height;
}

unsigned int gt_image_info_get_height(GtImageInfo *image_info)
{
  gt_assert(image_info);
  return image_info->height;
}

unsigned long gt_image_info_num_of_rec_maps(GtImageInfo *image_info)
{
  gt_assert(image_info);
  return gt_array_size(image_info->recmaps);
}

const GtRecMap* gt_image_info_get_rec_map(GtImageInfo *image_info,
                                          unsigned long i)
{
  gt_assert(image_info);
  return *(GtRecMap**) gt_array_get(image_info->recmaps, i);
}

int gt_image_info_unit_test(GtError *err)
{
  GtRecMap* rms[20];
  GtFeatureNode* features[20];
  GtImageInfo *ii;
  unsigned long i;
  GtStr *seqid;
  int had_err = 0;
  gt_assert(err);
  gt_error_check(err);

  seqid = gt_str_new_cstr("seqid");
  ii = gt_image_info_new();

  for (i=0;i<20;i++)
  {
    const GtRecMap* rm;
    unsigned long rbase;
    GtRange r;
    rbase = gt_rand_max(10);
    r.start = rbase;
    r.end = rbase + gt_rand_max(20);

    features[i] = (GtFeatureNode*) gt_feature_node_new(seqid, gt_ft_gene,
                                                       r.start, r.end,
                                                       GT_STRAND_FORWARD);
    rms[i] = gt_rec_map_new(gt_rand_max_double(100.0),
                            gt_rand_max_double(100.0),
                            gt_rand_max_double(100.0),
                            gt_rand_max_double(100.0),
                            features[i]);
    gt_image_info_add_rec_map(ii, rms[i]);
    gt_ensure(had_err, gt_image_info_num_of_rec_maps(ii) == i+1);
    gt_ensure(had_err, (rm = gt_image_info_get_rec_map(ii, i)) == rms[i]);
    gt_genome_node_delete((GtGenomeNode*) features[i]);
  }

  gt_image_info_delete(ii);
  gt_str_delete(seqid);

  return had_err;
}
