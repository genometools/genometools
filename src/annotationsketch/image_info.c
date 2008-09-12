/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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
#include "core/strand.h"
#include "core/unused_api.h"

struct GT_ImageInfo {
  GtArray* recmaps;
  unsigned int height;
};

GT_ImageInfo* gt_image_info_new()
{
  GT_ImageInfo *ii;
  ii = gt_calloc(1, sizeof (GT_ImageInfo));
  ii->recmaps = gt_array_new(sizeof (GT_RecMap*));
  assert(ii->recmaps);
  return ii;
}

void gt_image_info_delete(GT_ImageInfo *ii)
{
  unsigned long i;
  if (!ii) return;
  for (i=0;i<gt_image_info_num_of_recmaps(ii);i++)
  {
    GT_RecMap *rm = *(GT_RecMap**) gt_array_get(ii->recmaps, i);
    gt_recmap_delete(rm);
  }
  gt_array_delete(ii->recmaps);
  gt_free(ii);
}

void gt_image_info_add_recmap(GT_ImageInfo *ii, GT_RecMap *rm)
{
  assert(ii && rm);
  gt_array_add(ii->recmaps, rm);
}

void gt_image_info_set_height(GT_ImageInfo *ii, unsigned int height)
{
  assert(ii);
  ii->height = height;
}

unsigned int gt_image_info_get_height(GT_ImageInfo *ii)
{
  assert(ii);
  return ii->height;
}

unsigned long gt_image_info_num_of_recmaps(GT_ImageInfo *ii)
{
  assert(ii);
  return gt_array_size(ii->recmaps);
}

const GT_RecMap* gt_image_info_get_recmap(GT_ImageInfo *ii, unsigned long n)
{
  assert(ii);
  return *(GT_RecMap**) gt_array_get(ii->recmaps, n);
}

void gt_image_info_fill_recmap(GT_ImageInfo* ii, GT_RecMap* rm, unsigned long n)
{
  const GT_RecMap* tmp;
  assert(ii && rm);
  tmp = gt_image_info_get_recmap(ii, n);
  rm->nw_x = tmp->nw_x;
  rm->nw_y = tmp->nw_y;
  rm->se_x = tmp->se_x;
  rm->se_y = tmp->se_y;
  rm->gf = tmp->gf;
  rm->has_omitted_children = tmp->has_omitted_children;
}

int gt_image_info_unit_test(GT_Error *err)
{
  GT_RecMap* rms[20];
  GT_GenomeNode* gfs[20];
  GT_ImageInfo *ii;
  unsigned long i;
  GtStr *seqid;
  int had_err = 0;
  assert(err);
  gt_error_check(err);

  seqid = gt_str_new_cstr("seqid");
  ii = gt_image_info_new();

  for (i=0;i<20;i++)
  {
    const GT_RecMap* rm;
    unsigned long rbase;
    rbase = gt_rand_max(10);
    GT_Range r = { rbase, rbase + gt_rand_max(20)};
    gfs[i] = gt_genome_feature_new(seqid, gft_gene, r.start, r.end,
                                   GT_STRAND_FORWARD);
    rms[i] = gt_recmap_new(gt_rand_max_double(100.0),
                           gt_rand_max_double(100.0),
                           gt_rand_max_double(100.0),
                           gt_rand_max_double(100.0),
                           (GT_GenomeFeature*) /* XXX */ gfs[i]);
    gt_image_info_add_recmap(ii, rms[i]);
    ensure(had_err, gt_image_info_num_of_recmaps(ii) == i+1);
    ensure(had_err, (rm = gt_image_info_get_recmap(ii, i)) == rms[i]);
    ensure(had_err, rm->gf == rms[i]->gf);
    gt_genome_node_delete((GT_GenomeNode*) gfs[i]);
  }

  gt_image_info_delete(ii);
  gt_str_delete(seqid);

  return had_err;
}
