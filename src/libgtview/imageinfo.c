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

#include "libgtcore/array.h"
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/mathsupport.h"
#include "libgtcore/range.h"
#include "libgtcore/strand.h"
#include "libgtcore/unused.h"
#include "libgtext/feature_type_factory_builtin.h"
#include "libgtview/imageinfo.h"

struct ImageInfo {
  Array* recmaps;
  unsigned int height;
};

ImageInfo* image_info_new()
{
  ImageInfo *ii;
  ii = ma_calloc(1, sizeof (ImageInfo));
  ii->recmaps = array_new(sizeof (RecMap*));
  assert(ii->recmaps);
  return ii;
}

void image_info_delete(ImageInfo *ii)
{
  unsigned long i;
  if (!ii) return;
  for (i=0;i<image_info_num_of_recmaps(ii);i++)
  {
    RecMap *rm = *(RecMap**) array_get(ii->recmaps, i);
    recmap_delete(rm);
  }
  array_delete(ii->recmaps);
  ma_free(ii);
}

void image_info_add_recmap(ImageInfo *ii, RecMap *rm)
{
  assert(ii && rm);
  array_add(ii->recmaps, rm);
}

void image_info_set_height(ImageInfo *ii, unsigned int height)
{
  assert(ii);
  ii->height = height;
}

unsigned int image_info_get_height(ImageInfo *ii)
{
  assert(ii);
  return ii->height;
}

unsigned long image_info_num_of_recmaps(ImageInfo *ii)
{
  assert(ii);
  return array_size(ii->recmaps);
}

RecMap* image_info_get_recmap(ImageInfo *ii, unsigned long n)
{
  assert(ii);
  return *(RecMap**) array_get(ii->recmaps, n);
}

void image_info_get_recmap_ptr(ImageInfo *ii, RecMap *rm,  unsigned long n)
{
  RecMap* own_rm = *(RecMap**) array_get(ii->recmaps, n);
  assert(ii && rm);
  rm->nw_x = own_rm->nw_x;
  rm->nw_y = own_rm->nw_y;
  rm->se_x = own_rm->se_x;
  rm->se_y = own_rm->se_y;
  rm->gn = own_rm->gn;
}

int image_info_unit_test(Error *err)
{
  RecMap* rms[20];
  GenomeNode* gfs[20];
  FeatureTypeFactory *ftf;
  GenomeFeatureType *gft;
  ImageInfo *ii;
  unsigned long i;
  int had_err = 0;
  assert(err);
  error_check(err);

  ii = image_info_new();
  ftf = feature_type_factory_builtin_new();
  gft = feature_type_factory_create_gft(ftf, "gene");

  for (i=0;i<20;i++)
  {
    RecMap* rm;
    unsigned long rbase;
    rbase = rand_max(10);
    Range r = {rbase,rbase+rand_max(20)};
    gfs[i] = (GenomeNode*) genome_feature_new(gft, r, STRAND_FORWARD);
    rms[i] = recmap_create(rand_max_double(100.0),
                           rand_max_double(100.0),
                           rand_max_double(100.0),
                           rand_max_double(100.0),
                           gfs[i]);
    image_info_add_recmap(ii, rms[i]);
    ensure(had_err, image_info_num_of_recmaps(ii) == i+1);
    ensure(had_err, (rm = image_info_get_recmap(ii, i)) == rms[i]);
    ensure(had_err, rm->gn == rms[i]->gn);
    genome_node_delete((GenomeNode*) gfs[i]);
  }

  image_info_delete(ii);
  feature_type_factory_delete(ftf);
  return had_err;
}
