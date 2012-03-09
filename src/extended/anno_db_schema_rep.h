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

#ifndef ANNO_DB_SCHEMA_REP_H
#define ANNO_DB_SCHEMA_REP_H

#include "extended/anno_db_schema_api.h"
#include "extended/rdb_api.h"

typedef struct GtAnnoDBSchemaMembers GtAnnoDBSchemaMembers;
typedef struct GtAnnoDBSchemaClass GtAnnoDBSchemaClass;

struct GtAnnoDBSchema {
  const GtAnnoDBSchemaClass *c_class;
  GtAnnoDBSchemaMembers *members;
};

struct GtAnnoDBSchemaMembers {
  unsigned int reference_count;
};

typedef void            (*GtAnnoDBSchemaFreeFunc)(GtAnnoDBSchema*);
typedef GtFeatureIndex* (*GtAnnoDBSchemaBuildFunc)(GtAnnoDBSchema*, GtRDB*,
                                                   GtError*);

const GtAnnoDBSchemaClass* gt_anno_db_schema_class_new(size_t size,
                                               GtAnnoDBSchemaFreeFunc free_func,
                                               GtAnnoDBSchemaBuildFunc);

GtAnnoDBSchema* gt_anno_db_schema_create(const GtAnnoDBSchemaClass*);
void*           gt_anno_db_schema_cast(const GtAnnoDBSchemaClass*,
                                       GtAnnoDBSchema*);

#endif
