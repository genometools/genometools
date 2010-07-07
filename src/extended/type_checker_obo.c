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
#include "core/cstr_api.h"
#include "core/cstr_table.h"
#include "core/ma.h"
#include "extended/obo_parse_tree.h"
#include "extended/type_checker_obo.h"
#include "extended/type_checker_rep.h"

struct GtTypeCheckerOBO {
  const GtTypeChecker parent_instance;
  GtCstrTable *gt_feature_node_types;
};

#define gt_type_checker_obo_cast(FTF)\
        gt_type_checker_cast(gt_type_checker_obo_class(), FTF)

static void gt_type_checker_obo_free(GtTypeChecker *tc)
{
  GtTypeCheckerOBO *tco = gt_type_checker_obo_cast(tc);
  gt_cstr_table_delete(tco->gt_feature_node_types);
}

static bool gt_type_checker_obo_is_valid(GtTypeChecker *tc, const char *type)
{
  GtTypeCheckerOBO *tco;
  gt_assert(tc && type);
  tco = gt_type_checker_obo_cast(tc);
  return gt_cstr_table_get(tco->gt_feature_node_types, type) ? true : false;
}

const GtTypeCheckerClass* gt_type_checker_obo_class(void)
{
  static const GtTypeCheckerClass gt_type_checker_class =
    { sizeof (GtTypeCheckerOBO),
      gt_type_checker_obo_is_valid,
      gt_type_checker_obo_free };
  return &gt_type_checker_class;
}

static void add_feature_node_from_tree(GtTypeCheckerOBO *tco,
                                       GtOBOParseTree *obo_parse_tree,
                                       unsigned long stanza_num,
                                       const char *stanza_key)
{
  const char *value;
  gt_assert(tco && obo_parse_tree && stanza_key);
  value = gt_obo_parse_tree_get_stanza_value(obo_parse_tree, stanza_num,
                                          stanza_key);
  /* do not add values multiple times (possible for "name" values) */
  if (!gt_cstr_table_get(tco->gt_feature_node_types, value))
    gt_cstr_table_add(tco->gt_feature_node_types, value);
}

static int create_feature_nodes(GtTypeCheckerOBO *tco,
                                const char *obo_file_path, GtError *err)
{
  GtOBOParseTree *obo_parse_tree;
  unsigned long i;
  gt_error_check(err);
  gt_assert(tco && obo_file_path);
  if ((obo_parse_tree = gt_obo_parse_tree_new(obo_file_path, err))) {
    for (i = 0; i < gt_obo_parse_tree_num_of_stanzas(obo_parse_tree); i++) {
      if (!strcmp(gt_obo_parse_tree_get_stanza_type(obo_parse_tree, i),
                  "Term")) {
        const char *is_obsolete =
          gt_obo_parse_tree_get_stanza_value(obo_parse_tree, i, "is_obsolete");
        /* do not add obsolete types */
        if (!is_obsolete || strcmp(is_obsolete, "true")) {
          add_feature_node_from_tree(tco, obo_parse_tree, i, "id");
          add_feature_node_from_tree(tco, obo_parse_tree, i, "name");
        }
      }
    }
    gt_obo_parse_tree_delete(obo_parse_tree);
    return 0;
  }
  return -1;
}

GtTypeChecker* gt_type_checker_obo_new(const char *obo_file_path, GtError *err)
{
  GtTypeCheckerOBO *tco;
  GtTypeChecker *tc;
  gt_error_check(err);
  gt_assert(obo_file_path);
  tc = gt_type_checker_create(gt_type_checker_obo_class());
  tco = gt_type_checker_obo_cast(tc);
  tco->gt_feature_node_types = gt_cstr_table_new();
  if (create_feature_nodes(tco, obo_file_path, err)) {
    gt_type_checker_delete(tc);
    return NULL;
  }
  return tc;
}
