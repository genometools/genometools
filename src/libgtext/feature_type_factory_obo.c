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

#include "libgtcore/cstr.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtext/genome_feature_type_imp.h"
#include "libgtext/feature_type_factory_obo.h"
#include "libgtext/feature_type_factory_rep.h"
#include "libgtext/obo_parse_tree.h"

struct FeatureTypeFactoryOBO {
  const FeatureTypeFactory parent_instance;
  Hashtable *genome_feature_types;
};

#define feature_type_factory_obo_cast(FTF)\
        feature_type_factory_cast(feature_type_factory_obo_class(), FTF)

static void feature_type_factory_obo_free(FeatureTypeFactory *ftf)
{
  FeatureTypeFactoryOBO *ftfo = feature_type_factory_obo_cast(ftf);
  hashtable_delete(ftfo->genome_feature_types);
}

static GenomeFeatureType*
feature_type_factory_obo_create_gft(FeatureTypeFactory *ftf,
                                        const char *type)
{
  FeatureTypeFactoryOBO *ftfo;
  assert(ftf && type);
  ftfo = feature_type_factory_obo_cast(ftf);
  return hashtable_get(ftfo->genome_feature_types, type);
}

const FeatureTypeFactoryClass* feature_type_factory_obo_class(void)
{
  static const FeatureTypeFactoryClass feature_type_factory_class =
    { sizeof (FeatureTypeFactoryOBO),
      feature_type_factory_obo_create_gft,
      feature_type_factory_obo_free };
  return &feature_type_factory_class;
}

static void add_genome_feature_from_tree(FeatureTypeFactoryOBO *ftfo,
                                         OBOParseTree *obo_parse_tree,
                                         unsigned long stanza_num,
                                         const char *stanza_key)
{
  GenomeFeatureType *type;
  char *value;
  assert(ftfo && obo_parse_tree && stanza_key);
  value = cstr_dup(obo_parse_tree_get_stanza_value(obo_parse_tree, stanza_num,
                                                   stanza_key));
  type = genome_feature_type_construct((FeatureTypeFactory*) ftfo, value);
  hashtable_add(ftfo->genome_feature_types, value, type);
}

static int create_genome_features(FeatureTypeFactoryOBO *ftfo,
                                  const char *obo_file_path, Error *err)
{
  OBOParseTree *obo_parse_tree;
  unsigned long i;
  error_check(err);
  assert(ftfo && obo_file_path);
  if ((obo_parse_tree = obo_parse_tree_new(obo_file_path, err))) {
    for (i = 0; i < obo_parse_tree_num_of_stanzas(obo_parse_tree); i++) {
     add_genome_feature_from_tree(ftfo, obo_parse_tree, i, "id");
     add_genome_feature_from_tree(ftfo, obo_parse_tree, i, "name");
    }
    obo_parse_tree_delete(obo_parse_tree);
    return 0;
  }
  return -1;
}

FeatureTypeFactory* feature_type_factory_obo_new(const char *obo_file_path,
                                                 Error *err)
{
  FeatureTypeFactoryOBO *ftfo;
  FeatureTypeFactory *ftf;
  error_check(err);
  assert(obo_file_path);
  ftf = feature_type_factory_create(feature_type_factory_obo_class());
  ftfo = feature_type_factory_obo_cast(ftf);
  ftfo->genome_feature_types = hashtable_new(HASH_STRING, ma_free_func,
                                             (FreeFunc)
                                             genome_feature_type_delete);
  if (create_genome_features(ftfo, obo_file_path, err)) {
    feature_type_factory_delete(ftf);
    return NULL;
  }
  return ftf;
}
