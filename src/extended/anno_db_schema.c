/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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
#include "core/ma.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/anno_db_schema_api.h"
#include "extended/anno_db_schema_rep.h"

struct GtAnnoDBSchemaClass {
  size_t size;
  GtAnnoDBSchemaFreeFunc free_func;
  GtAnnoDBSchemaBuildFunc build_func;
};

const GtAnnoDBSchemaClass* gt_anno_db_schema_class_new(size_t size,
                            GtAnnoDBSchemaFreeFunc free_func,
                            GtAnnoDBSchemaBuildFunc build_func)
{
  GtAnnoDBSchemaClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free_func = free_func;
  c_class->build_func = build_func;
  return c_class;
}

GtAnnoDBSchema* gt_anno_db_schema_create(const GtAnnoDBSchemaClass *dbc)
{
  GtAnnoDBSchema *s;
  gt_assert(dbc && dbc->size);
  s = gt_calloc(1, dbc->size);
  s->c_class = dbc;
  s->members = gt_calloc(1, sizeof (GtAnnoDBSchemaMembers));
  return s;
}

GtAnnoDBSchema* gt_anno_db_schema_ref(GtAnnoDBSchema *s)
{
  gt_assert(s);
  s->members->reference_count++;
  return s;
}

void gt_anno_db_schema_delete(GtAnnoDBSchema *s)
{
  if (!s) return;
  if (s->members->reference_count) {
    s->members->reference_count--;
    return;
  }
  gt_assert(s->c_class);
  if (s->c_class->free_func) s->c_class->free_func(s);
  gt_free(s->members);
  gt_free(s);
}

void* gt_anno_db_schema_cast(GT_UNUSED const GtAnnoDBSchemaClass *sc,
                             GtAnnoDBSchema *s)
{
  gt_assert(sc && s && s->c_class == sc);
  return s;
}

GtFeatureIndex* gt_anno_db_schema_get_feature_index(GtAnnoDBSchema *s,
                                                    GtRDB *rdb, GtError *err)
{
  gt_assert(s && s->c_class);
  if (s->c_class->build_func)
    return s->c_class->build_func(s, rdb, err);
  return NULL;
}
