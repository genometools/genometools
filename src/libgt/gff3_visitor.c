/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "fptr.h"
#include "genome_node.h"
#include "genome_visitor_rep.h"
#include "gff3_output.h"
#include "gff3_parser.h"
#include "gff3_visitor.h"
#include "hashtable.h"
#include "str.h"
#include "undef.h"
#include "xansi.h"

struct Gff3_visitor {
  const GenomeVisitor parent_instance;
  unsigned int version_string_shown;
  unsigned long *id_counter;
  Hashtable *genome_feature_to_id_array,
            *genome_feature_to_unique_id_str;
  FILE *outfp;
};

typedef struct {
  Hashtable *genome_feature_to_id_array;
  const char *id;
} Add_id_info;

#define gff3_visitor_cast(GV)\
        genome_visitor_cast(gff3_visitor_class(), GV)

static void gff3_version_string(GenomeVisitor *gv)
{
  Gff3_visitor *gff3_visitor = gff3_visitor_cast(gv);
  assert(gff3_visitor);
  if (!gff3_visitor->version_string_shown) {
    fprintf(gff3_visitor->outfp, "%s\n", GFF_VERSION_STRING);
    gff3_visitor->version_string_shown = 1;
  }
}

static void gff3_visitor_free(GenomeVisitor *gv)
{
  Gff3_visitor *gff3_visitor = gff3_visitor_cast(gv);
  assert(gff3_visitor);
  free(gff3_visitor->id_counter);
  hashtable_free(gff3_visitor->genome_feature_to_id_array);
  hashtable_free(gff3_visitor->genome_feature_to_unique_id_str);
}

static int gff3_visitor_comment(GenomeVisitor *gv, Comment *c,
                                /*@unused@*/ Log *l, Error *err)
{
  Gff3_visitor *gff3_visitor;
  error_check(err);
  gff3_visitor = gff3_visitor_cast(gv);
  assert(gv && c);
  gff3_version_string(gv);
  fprintf(gff3_visitor->outfp, "#%s\n", comment_get_comment(c));
  return 0;
}

static int add_id(GenomeNode *gn, void *data, Error *err)
{
  Add_id_info *info = (Add_id_info*) data;
  Array *parent_features = NULL;
  error_check(err);
  assert(gn && info && info->genome_feature_to_id_array && info->id);
  parent_features = hashtable_get(info->genome_feature_to_id_array, gn);
  if (!parent_features) {
    parent_features = array_new(sizeof(char*));
    hashtable_add(info->genome_feature_to_id_array, gn, parent_features);
  }
  array_add(parent_features, info->id);
  return 0;
}

static int gff3_show_genome_feature(GenomeNode *gn, void *data, Error *err)
{
  unsigned int part_shown = 0;
  Gff3_visitor *gff3_visitor = (Gff3_visitor*) data;
  Genome_feature *gf = (Genome_feature*) gn;
  Array *parent_features = NULL;
  unsigned long i;
  Str *id;

  error_check(err);
  assert(gn && gf && gff3_visitor);

  /* output leading part */
  gff3_output_leading(gf, gff3_visitor->outfp);

  /* show unique id part of attributes */
  if ((id = hashtable_get(gff3_visitor->genome_feature_to_unique_id_str, gn))) {
    fprintf(gff3_visitor->outfp, "%s=%s", ID_STRING, str_get(id));
    part_shown = 1;
  }

  /* show parent part of attributes */
  parent_features = hashtable_get(gff3_visitor->genome_feature_to_id_array,
                                  gn);
  if (array_size(parent_features)) {
    if (part_shown)
      xfputc(';', gff3_visitor->outfp);
    fprintf(gff3_visitor->outfp, "%s=", PARENT_STRING);
    for (i = 0; i < array_size(parent_features); i++) {
      if (i)
        xfputc(',', gff3_visitor->outfp);
      fprintf(gff3_visitor->outfp, "%s",
              *(char**) array_get(parent_features, i));
    }
    part_shown = 1;
  }

  /* XXX: show missing part of attributes */
  /* if (part_shown) putc(';', gff3_visitor->outfp); */

  xfputc('\n', gff3_visitor->outfp);

  return 0;
}

static int store_ids(GenomeNode *gn, void *data, Error *err)
{
  Gff3_visitor *gff3_visitor = (Gff3_visitor*) data;
  Genome_feature *gf = (Genome_feature*) gn;
  GenomeFeatureType type;
  Add_id_info add_id_info;
  int has_err = 0;
  Str *id;

  error_check(err);
  assert(gn && gf && gff3_visitor);
  type = genome_feature_get_type(gf);

  if (genome_node_has_children(gn)) {
    /* increase id counter */
    gff3_visitor->id_counter[type]++;

    /* build id string */
    id = str_new_cstr(genome_feature_type_get_cstr(type));
    str_append_ulong(id, gff3_visitor->id_counter[type]);

    /* store (unique) id */
    hashtable_add(gff3_visitor->genome_feature_to_unique_id_str, gn, id);

    /* for each child -> store the parent feature in the hash table */
    add_id_info.genome_feature_to_id_array =
      gff3_visitor->genome_feature_to_id_array,
    add_id_info.id = str_get(id);
    has_err = genome_node_traverse_direct_children(gn, &add_id_info, add_id,
                                                   err);
  }
  return has_err;
}

static int gff3_visitor_genome_feature(GenomeVisitor *gv, Genome_feature *gf,
                                       /*@unused@*/ Log *l, Error *err)
{
  Gff3_visitor *gff3_visitor;
  int has_err;
  error_check(err);
  gff3_visitor = gff3_visitor_cast(gv);

  gff3_version_string(gv);

  has_err = genome_node_traverse_children((GenomeNode*) gf, gff3_visitor,
                                          store_ids, true, err);
  if (!has_err) {
    has_err = genome_node_traverse_children((GenomeNode*) gf, gff3_visitor,
                                            gff3_show_genome_feature, true,
                                            err);
  }

  /* clear hashtable */
  /* XXX */
  hashtable_reset(gff3_visitor->genome_feature_to_id_array);
  hashtable_reset(gff3_visitor->genome_feature_to_unique_id_str);

  /* terminator */
  xfputs("###\n", gff3_visitor->outfp);

  return has_err;
}

static int gff3_visitor_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                        /*@unused@*/ Log *l, Error *err)
{
  Gff3_visitor *gff3_visitor;
  error_check(err);
  gff3_visitor = gff3_visitor_cast(gv);
  assert(gv && sr);
  /* a sequence region has no children */
  assert(!genome_node_has_children((GenomeNode*) sr));

  gff3_version_string(gv);
  fprintf(gff3_visitor->outfp,
          "%s   %s %lu %lu\n",
          GFF_SEQUENCE_REGION,
          str_get(genome_node_get_seqid((GenomeNode*) sr)),
          genome_node_get_start((GenomeNode*) sr),
          genome_node_get_end((GenomeNode*) sr));
  return 0;
}

const GenomeVisitorClass* gff3_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof(Gff3_visitor),
                                          gff3_visitor_free,
                                          gff3_visitor_comment,
                                          gff3_visitor_genome_feature,
                                          gff3_visitor_sequence_region,
                                          NULL };
  return &gvc;
}

GenomeVisitor* gff3_visitor_new(FILE *outfp)
{
  GenomeVisitor *gv = genome_visitor_create(gff3_visitor_class());
  Gff3_visitor *gff3_visitor = gff3_visitor_cast(gv);
  gff3_visitor->version_string_shown = 0;
  gff3_visitor->id_counter = xcalloc(genome_feature_type_num_of_features(),
                                     sizeof(unsigned long));
  gff3_visitor->genome_feature_to_id_array = hashtable_new(HASH_DIRECT, NULL,
                                                           (Free) array_free);
  gff3_visitor->genome_feature_to_unique_id_str = hashtable_new(HASH_DIRECT,
                                                                NULL,
                                                                (Free)
                                                                str_free);
  gff3_visitor->outfp = outfp;
  return gv;
}
