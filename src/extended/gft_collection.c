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

#include "core/cstr.h"
#include "core/cstr_table.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/gft_collection.h"

struct GFTCollection {
  GT_CstrTable *genome_feature_types;
};

GFTCollection* gft_collection_new(void)
{
  GFTCollection *gftc = gt_malloc(sizeof (GFTCollection));
  gftc->genome_feature_types = gt_cstr_table_new();
  return gftc;
}

void gft_collection_delete(GFTCollection *gftc)
{
  if (!gftc) return;
  gt_cstr_table_delete(gftc->genome_feature_types);
  gt_free(gftc);
}

void gft_collection_add(GFTCollection *gftc, const char *type)
{
  assert(gftc && type);
  gt_cstr_table_add(gftc->genome_feature_types, type);
}

const char* gft_collection_get(GFTCollection *gftc, const char *type)
{
  assert(gftc && type);
  return gt_cstr_table_get(gftc->genome_feature_types, type);
}

GT_StrArray* gft_collection_get_types(const GFTCollection *gftc)
{
  assert(gftc);
  return gt_cstr_table_get_all(gftc->genome_feature_types);
}
