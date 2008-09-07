/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include <stdlib.h>
#include <string.h>
#include "core/fasta.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/unused.h"
#include "core/string_distri.h"
#include "extended/genome_node.h"
#include "extended/genome_visitor_rep.h"
#include "extended/gff3_output.h"
#include "extended/gff3_parser.h"
#include "extended/gff3_visitor.h"

struct GFF3Visitor {
  const GenomeVisitor parent_instance;
  bool version_string_shown,
       fasta_directive_shown;
  StringDistri *id_counter;
  Hashmap *gt_genome_feature_to_id_array,
            *gt_genome_feature_to_unique_id_str;
  unsigned long fasta_width;
  GenFile *outfp;
};

typedef struct {
  Hashmap *gt_genome_feature_to_id_array;
  const char *id;
} Add_id_info;

typedef struct {
  bool *attribute_shown;
  GenFile *outfp;
} ShowAttributeInfo;

#define gff3_visitor_cast(GV)\
        genome_visitor_cast(gff3_visitor_class(), GV)

static void gff3_version_string(GenomeVisitor *gv)
{
  GFF3Visitor *gff3_visitor = gff3_visitor_cast(gv);
  assert(gff3_visitor);
  if (!gff3_visitor->version_string_shown) {
    genfile_xprintf(gff3_visitor->outfp, "%s   %u\n", GFF_VERSION_PREFIX,
                    GFF_VERSION);
    gff3_visitor->version_string_shown = true;
  }
}

static void gff3_visitor_free(GenomeVisitor *gv)
{
  GFF3Visitor *gff3_visitor = gff3_visitor_cast(gv);
  assert(gff3_visitor);
  string_distri_delete(gff3_visitor->id_counter);
  hashmap_delete(gff3_visitor->gt_genome_feature_to_id_array);
  hashmap_delete(gff3_visitor->gt_genome_feature_to_unique_id_str);
}

static int gff3_visitor_comment(GenomeVisitor *gv, GT_Comment *c,
                                UNUSED GT_Error *err)
{
  GFF3Visitor *gff3_visitor;
  gt_error_check(err);
  gff3_visitor = gff3_visitor_cast(gv);
  assert(gv && c);
  gff3_version_string(gv);
  genfile_xprintf(gff3_visitor->outfp, "#%s\n", gt_comment_get_comment(c));
  return 0;
}

static int add_id(GT_GenomeNode *gn, void *data, UNUSED GT_Error *err)
{
  Add_id_info *info = (Add_id_info*) data;
  GT_Array *parent_features = NULL;
  gt_error_check(err);
  assert(gn && info && info->gt_genome_feature_to_id_array && info->id);
  parent_features = hashmap_get(info->gt_genome_feature_to_id_array, gn);
  if (!parent_features) {
    parent_features = gt_array_new(sizeof (char*));
    hashmap_add(info->gt_genome_feature_to_id_array, gn, parent_features);
  }
  gt_array_add(parent_features, info->id);
  return 0;
}

static void show_attribute(const char *attr_name, const char *attr_value,
                           void *data)
{
  ShowAttributeInfo *info = (ShowAttributeInfo*) data;
  assert(attr_name && attr_value && info);
  if (strcmp(attr_name, ID_STRING) && strcmp(attr_name, PARENT_STRING)) {
    if (*info->attribute_shown)
      genfile_xfputc(';', info->outfp);
    else
      *info->attribute_shown = true;
    genfile_xprintf(info->outfp, "%s=%s", attr_name, attr_value);
  }
}

static int gff3_show_genome_feature(GT_GenomeNode *gn, void *data,
                                    UNUSED GT_Error *err)
{
  bool part_shown = false;
  GFF3Visitor *gff3_visitor = (GFF3Visitor*) data;
  GT_GenomeFeature *gf = (GT_GenomeFeature*) gn;
  GT_Array *parent_features = NULL;
  ShowAttributeInfo info;
  unsigned long i;
  GT_Str *id;

  gt_error_check(err);
  assert(gn && gf && gff3_visitor);

  /* output leading part */
  gff3_output_leading(gf, gff3_visitor->outfp);

  /* show unique id part of attributes */
  if ((id = hashmap_get(gff3_visitor->gt_genome_feature_to_unique_id_str, gn))) {
    genfile_xprintf(gff3_visitor->outfp, "%s=%s", ID_STRING, gt_str_get(id));
    part_shown = true;
  }

  /* show parent part of attributes */
  parent_features = hashmap_get(gff3_visitor->gt_genome_feature_to_id_array,
                                  gn);
  if (gt_array_size(parent_features)) {
    if (part_shown)
      genfile_xfputc(';', gff3_visitor->outfp);
    genfile_xprintf(gff3_visitor->outfp, "%s=", PARENT_STRING);
    for (i = 0; i < gt_array_size(parent_features); i++) {
      if (i)
        genfile_xfputc(',', gff3_visitor->outfp);
      genfile_xprintf(gff3_visitor->outfp, "%s",
                      *(char**) gt_array_get(parent_features, i));
    }
    part_shown = true;
  }

  /* show missing part of attributes */
  info.attribute_shown = &part_shown;
  info.outfp = gff3_visitor->outfp;
  gt_genome_feature_foreach_attribute(gf, show_attribute, &info);

  /* show dot if no attributes have been shown */
  if (!part_shown)
    genfile_xfputc('.', gff3_visitor->outfp);

  /* show terminal newline */
  genfile_xfputc('\n', gff3_visitor->outfp);

  return 0;
}

static GT_Str* create_unique_id(GFF3Visitor *gff3_visitor, GT_GenomeFeature *gf)
{
  GT_GenomeFeatureType *type;
  GT_Str *id;
  assert(gff3_visitor && gf);
  type = gt_genome_feature_get_type(gf);

  /* increase id counter */
  string_distri_add(gff3_visitor->id_counter,
                   gt_genome_feature_type_get_cstr(type));

  /* build id string */
  id = gt_str_new_cstr(gt_genome_feature_type_get_cstr(type));
  gt_str_append_ulong(id, string_distri_get(gff3_visitor->id_counter,
                                        gt_genome_feature_type_get_cstr(type)));
  /* store (unique) id */
  hashmap_add(gff3_visitor->gt_genome_feature_to_unique_id_str, gf, id);

  return id;
}

static int store_ids(GT_GenomeNode *gn, void *data, GT_Error *err)
{
  GFF3Visitor *gff3_visitor = (GFF3Visitor*) data;
  GT_GenomeFeature *gf = (GT_GenomeFeature*) gn;
  Add_id_info add_id_info;
  int had_err = 0;
  GT_Str *id;

  gt_error_check(err);
  assert(gn && gf && gff3_visitor);

  if (gt_genome_node_has_children(gn) || gt_genome_feature_is_multi(gf)) {
    if (gt_genome_feature_is_multi(gf)) {
      id = hashmap_get(gff3_visitor->gt_genome_feature_to_unique_id_str,
                       gt_genome_feature_get_multi_representative(gf));
      if (!id) { /* the representative does not have its own id */
        id = create_unique_id(gff3_visitor,
                              gt_genome_feature_get_multi_representative(gf));
      }
      if (gt_genome_feature_get_multi_representative(gf) != gf) {
        hashmap_add(gff3_visitor->gt_genome_feature_to_unique_id_str, gf,
                    gt_str_ref(id));
      }
    }
    else
      id = create_unique_id(gff3_visitor, gf);

    /* for each child -> store the parent feature in the hash map */
    add_id_info.gt_genome_feature_to_id_array =
      gff3_visitor->gt_genome_feature_to_id_array,
    add_id_info.id = gt_str_get(id);
    had_err = gt_genome_node_traverse_direct_children(gn, &add_id_info, add_id,
                                                   err);
  }
  return had_err;
}

static int gff3_visitor_genome_feature(GenomeVisitor *gv, GT_GenomeFeature *gf,
                                       GT_Error *err)
{
  GFF3Visitor *gff3_visitor;
  int had_err;
  gt_error_check(err);
  gff3_visitor = gff3_visitor_cast(gv);

  gff3_version_string(gv);

  had_err = gt_genome_node_traverse_children((GT_GenomeNode*) gf, gff3_visitor,
                                          store_ids, true, err);
  if (!had_err) {
    if (gt_genome_node_is_tree((GT_GenomeNode*) gf)) {
      had_err = gt_genome_node_traverse_children((GT_GenomeNode*) gf, gff3_visitor,
                                              gff3_show_genome_feature, true,
                                              err);
    }
    else {
      /* got a DAG -> traverse bin breadth first fashion to make sure that the
         'Parent' attributes are shown in correct order */
      had_err = gt_genome_node_traverse_children_breadth((GT_GenomeNode*) gf,
                                                      gff3_visitor,
                                                      gff3_show_genome_feature,
                                                      true, err);
    }
  }

  /* reset hashmaps */
  hashmap_reset(gff3_visitor->gt_genome_feature_to_id_array);
  hashmap_reset(gff3_visitor->gt_genome_feature_to_unique_id_str);

  /* show terminator, if the feature has children (otherwise it is clear that
     the feature is complete, because no ID attribute has been shown) */
  if (gt_genome_node_has_children((GT_GenomeNode*) gf))
    genfile_xprintf(gff3_visitor->outfp, "%s\n", GFF_TERMINATOR);

  return had_err;
}

static int gff3_visitor_sequence_region(GenomeVisitor *gv, GT_SequenceRegion *sr,
                                        UNUSED GT_Error *err)
{
  GFF3Visitor *gff3_visitor;
  gt_error_check(err);
  gff3_visitor = gff3_visitor_cast(gv);
  assert(gv && sr);
  /* a sequence region has no children */
  assert(!gt_genome_node_has_children((GT_GenomeNode*) sr));

  gff3_version_string(gv);
  genfile_xprintf(gff3_visitor->outfp, "%s   %s %lu %lu\n", GFF_SEQUENCE_REGION,
                  gt_str_get(gt_genome_node_get_seqid((GT_GenomeNode*) sr)),
                  gt_genome_node_get_start((GT_GenomeNode*) sr),
                  gt_genome_node_get_end((GT_GenomeNode*) sr));
  return 0;
}

static int gff3_visitor_sequence_node(GenomeVisitor *gv, GT_SequenceNode *sn,
                                      UNUSED GT_Error *err)
{
  GFF3Visitor *gff3_visitor;
  gt_error_check(err);
  gff3_visitor = gff3_visitor_cast(gv);
  assert(gv && sn);
  if (!gff3_visitor->fasta_directive_shown) {
    genfile_xprintf(gff3_visitor->outfp, "%s\n", GFF_FASTA_DIRECTIVE);
    gff3_visitor->fasta_directive_shown = true;
  }
  fasta_show_entry_generic(gt_sequence_node_get_description(sn),
                           gt_sequence_node_get_sequence(sn),
                           gt_sequence_node_get_sequence_length(sn),
                           gff3_visitor->fasta_width, gff3_visitor->outfp);
  return 0;
}

const GenomeVisitorClass* gff3_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (GFF3Visitor),
                                          gff3_visitor_free,
                                          gff3_visitor_comment,
                                          gff3_visitor_genome_feature,
                                          gff3_visitor_sequence_region,
                                          gff3_visitor_sequence_node };
  return &gvc;
}

GenomeVisitor* gff3_visitor_new(GenFile *outfp)
{
  GenomeVisitor *gv = genome_visitor_create(gff3_visitor_class());
  GFF3Visitor *gff3_visitor = gff3_visitor_cast(gv);
  gff3_visitor->version_string_shown = false;
  gff3_visitor->fasta_directive_shown = false;
  gff3_visitor->id_counter = string_distri_new();
  gff3_visitor->gt_genome_feature_to_id_array = hashmap_new(
    HASH_DIRECT, NULL, (FreeFunc) gt_array_delete);
  gff3_visitor->gt_genome_feature_to_unique_id_str = hashmap_new(
    HASH_DIRECT, NULL, (FreeFunc) gt_str_delete);
  gff3_visitor->fasta_width = 0;
  gff3_visitor->outfp = outfp;
  return gv;
}

void gff3_visitor_set_fasta_width(GenomeVisitor *gv, unsigned long fasta_width)
{
  GFF3Visitor *gff3_visitor = gff3_visitor_cast(gv);
  assert(gv);
  gff3_visitor->fasta_width = fasta_width;
}
