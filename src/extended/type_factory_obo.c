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

#include <string.h>
#include "core/cstr.h"
#include "core/cstr_table.h"
#include "core/ma.h"
#include "extended/feature_type_imp.h"
#include "extended/obo_parse_tree.h"
#include "extended/type_factory_obo.h"
#include "extended/type_factory_rep.h"

struct GT_TypeFactoryOBO {
  const GT_TypeFactory parent_instance;
  CstrTable *gt_genome_feature_types;
};

#define gt_type_factory_obo_cast(FTF)\
        gt_type_factory_cast(gt_type_factory_obo_class(), FTF)

static void gt_type_factory_obo_free(GT_TypeFactory *ftf)
{
  GT_TypeFactoryOBO *ftfo = gt_type_factory_obo_cast(ftf);
  cstr_table_delete(ftfo->gt_genome_feature_types);
}

static GT_FeatureType* gt_type_factory_obo_create_gft(GT_TypeFactory *ftf,
                                                      const char *type)
{
  GT_TypeFactoryOBO *ftfo;
  GT_FeatureType *gft;
  assert(ftf && type);
  ftfo = gt_type_factory_obo_cast(ftf);
  if (!(gft = gft_collection_get(ftf->used_types, type))) {
    if (cstr_table_get(ftfo->gt_genome_feature_types, type)) {
      gft = gt_feature_type_construct(ftf, type);
      gft_collection_add(ftf->used_types, type, gft);
    }
  }
  return gft;
}

const GT_TypeFactoryClass* gt_type_factory_obo_class(void)
{
  static const GT_TypeFactoryClass gt_type_factory_class =
    { sizeof (GT_TypeFactoryOBO),
      gt_type_factory_obo_create_gft,
      gt_type_factory_obo_free };
  return &gt_type_factory_class;
}

static void add_gt_genome_feature_from_tree(GT_TypeFactoryOBO *ftfo,
                                         OBOParseTree *obo_parse_tree,
                                         unsigned long stanza_num,
                                         const char *stanza_key)
{
  const char *value;
  assert(ftfo && obo_parse_tree && stanza_key);
  value = obo_parse_tree_get_stanza_value(obo_parse_tree, stanza_num,
                                          stanza_key);
  /* do not add values multiple times (possible for "name" values) */
  if (!cstr_table_get(ftfo->gt_genome_feature_types, value))
    cstr_table_add(ftfo->gt_genome_feature_types, value);
}

static int create_genome_features(GT_TypeFactoryOBO *ftfo,
                                  const char *obo_file_path, GT_Error *err)
{
  OBOParseTree *obo_parse_tree;
  unsigned long i;
  gt_error_check(err);
  assert(ftfo && obo_file_path);
  if ((obo_parse_tree = obo_parse_tree_new(obo_file_path, err))) {
    for (i = 0; i < obo_parse_tree_num_of_stanzas(obo_parse_tree); i++) {
      if (!strcmp(obo_parse_tree_get_stanza_type(obo_parse_tree, i), "Term")) {
        const char *is_obsolete =
          obo_parse_tree_get_stanza_value(obo_parse_tree, i, "is_obsolete");
        /* do not add obsolete types */
        if (!is_obsolete || strcmp(is_obsolete, "true")) {
          add_gt_genome_feature_from_tree(ftfo, obo_parse_tree, i, "id");
          add_gt_genome_feature_from_tree(ftfo, obo_parse_tree, i, "name");
        }
      }
    }
    obo_parse_tree_delete(obo_parse_tree);
    return 0;
  }
  return -1;
}

GT_TypeFactory* gt_type_factory_obo_new(const char *obo_file_path,
                                                 GT_Error *err)
{
  GT_TypeFactoryOBO *ftfo;
  GT_TypeFactory *ftf;
  gt_error_check(err);
  assert(obo_file_path);
  ftf = gt_type_factory_create(gt_type_factory_obo_class());
  ftfo = gt_type_factory_obo_cast(ftf);
  ftfo->gt_genome_feature_types = cstr_table_new();
  if (create_genome_features(ftfo, obo_file_path, err)) {
    gt_type_factory_delete(ftf);
    return NULL;
  }
  return ftf;
}
