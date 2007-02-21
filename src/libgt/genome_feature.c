/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include "genome_feature.h"
#include "genome_feature_type.h"
#include "genome_node_rep.h"
#include "hashtable.h"
#include "range.h"
#include "strand.h"
#include "undef.h"
#include "xansi.h"

struct GenomeFeature
{
  const GenomeNode parent_instance;
  Str *seqid,
      *source;
  GenomeFeatureType type;
  Range range;
  double score;
  Strand strand;
  Phase phase;
  Hashtable *attributes; /* stores additional the attributes besides 'Parent'
                            and 'ID'; created on demand */
};

#define genome_feature_cast(GN)\
        genome_node_cast(genome_feature_class(), GN)

static void genome_feature_free(GenomeNode *gn)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  assert(gf);
  str_delete(gf->seqid);
  str_delete(gf->source);
  hashtable_delete(gf->attributes);
}

static Str* genome_feature_get_seqid(GenomeNode *gn)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  return gf->seqid;
}

static Range genome_feature_get_range(GenomeNode *gn)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  return gf->range;
}

static void genome_feature_set_seqid(GenomeNode *gn, Str *seqid)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  assert(gf && seqid && !gf->seqid);
  gf->seqid = str_ref(seqid);
}

static void genome_feature_set_source(GenomeNode *gn, Str *source)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  assert(gf && source && !gf->source);
  gf->source = str_ref(source);
}

static void genome_feature_set_phase(GenomeNode *gn, Phase phase)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  assert(gf && gf->phase == PHASE_UNDEFINED);
  gf->phase = phase;
}

static int genome_feature_accept(GenomeNode *gn, GenomeVisitor *gv, Env *env)
{
  GenomeFeature *gf;
  env_error_check(env);
  gf = genome_feature_cast(gn);
  return genome_visitor_visit_genome_feature(gv, gf, env);
}

const GenomeNodeClass* genome_feature_class()
{
  static const GenomeNodeClass gnc = { sizeof (GenomeFeature),
                                       genome_feature_free,
                                       genome_feature_get_seqid,
                                       genome_feature_get_seqid,
                                       genome_feature_get_range,
                                       NULL,
                                       genome_feature_set_seqid,
                                       genome_feature_set_source,
                                       genome_feature_set_phase,
                                       genome_feature_accept };
  return &gnc;
}

GenomeNode* genome_feature_new(GenomeFeatureType type,
                                Range range,
                                Strand strand,
                                const char *filename,
                                unsigned long line_number)
{
  GenomeNode *gn = genome_node_create(genome_feature_class(), filename,
                                       line_number);
  GenomeFeature *gf = genome_feature_cast(gn);
  assert(range.start <= range.end);
  gf->seqid      = NULL;
  gf->source     = NULL;
  gf->type       = type;
  gf->score      = UNDEFDOUBLE;
  gf->range      = range;
  gf->strand     = strand;
  gf->phase      = PHASE_UNDEFINED;
  gf->attributes = NULL;
  return gn;
}

const char* genome_feature_get_source(GenomeFeature *gf)
{
  assert(gf);
  return gf->source ? str_get(gf->source) : ".";
}

GenomeFeatureType genome_feature_get_type(GenomeFeature *gf)
{
  assert(gf);
  return gf->type;
}

double genome_feature_get_score(GenomeFeature *gf)
{
  assert(gf);
  return gf->score;
}

Strand genome_feature_get_strand(GenomeFeature *gf)
{
  assert(gf);
  return gf->strand;
}

Phase genome_feature_get_phase(GenomeFeature *gf)
{
  assert(gf);
  return gf->phase;
}

static int save_exon(GenomeNode *gn, void *data, Env *env)
{
  GenomeFeature *gf;
  Array *exon_features = (Array*) data;
  env_error_check(env);
  gf = (GenomeFeature*) gn;
  assert(gf && exon_features);
  if (genome_feature_get_type(gf) == gft_exon) {
    array_add(exon_features, gf);
  }
  return 0;
}

void genome_feature_get_exons(GenomeFeature *gf, Array *exon_features)
{
  int has_err;
  assert(gf && exon_features && !array_size(exon_features));
  has_err = genome_node_traverse_children((GenomeNode*) gf, exon_features,
                                          save_exon, false, NULL);
  assert(!has_err); /* cannot happen, because save_exon() is sane */
}

void genome_feature_set_end(GenomeFeature *gf, unsigned long end)
{
  assert(gf && gf->range.start <= end);
  gf->range.end = end;
}

void genome_feature_set_score(GenomeFeature *gf, double score)
{
  assert(gf);
  gf->score = score;
}

void genome_feature_add_attribute(GenomeFeature *gf, const char *attr_name,
                                  const char *attr_value)
{
  assert(gf && attr_name && attr_value);
  if (!gf->attributes)
    gf->attributes = hashtable_new(HASH_DIRECT, free, free);
  hashtable_add(gf->attributes, xstrdup(attr_name), xstrdup(attr_value));
}

bool genome_feature_has_attribute(const GenomeFeature *gf)
{
  assert(gf);
  if (gf->attributes)
    return true;
  return false;
}

int genome_feature_foreach_attribute(GenomeFeature *gf,
                                     AttributeIterFunc iterfunc, void *data,
                                     Env *env)
{
  env_error_check(env);
  assert(gf && iterfunc);
  assert(genome_feature_has_attribute(gf));
  return hashtable_foreach(gf->attributes, (Hashiteratorfunc) iterfunc, data,
                           env);
}
