/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include <libgt/cstr.h>
#include <libgt/genome_feature.h>
#include <libgt/genome_feature_type.h>
#include <libgt/genome_node_rep.h>
#include <libgt/hashtable.h>
#include <libgt/range.h>
#include <libgt/strand.h>
#include <libgt/undef.h>
#include <libgt/xansi.h>

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
  TranscriptFeatureType transcripttype;
};

typedef struct {
  Array *exon_features,
        *cds_features;
} SaveExonAndCDSInfo;

#define genome_feature_cast(GN)\
        genome_node_cast(genome_feature_class(), GN)

static void genome_feature_free(GenomeNode *gn, Env *env)
{
  GenomeFeature *gf = genome_feature_cast(gn);
  assert(gf);
  str_delete(gf->seqid, env);
  str_delete(gf->source, env);
  hashtable_delete(gf->attributes, env);
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
                                unsigned long line_number, Env *env)
{
  GenomeNode *gn = genome_node_create(genome_feature_class(), filename,
                                      line_number, env);
  GenomeFeature *gf = genome_feature_cast(gn);
  assert(range.start <= range.end);
  gf->seqid          = NULL;
  gf->source         = NULL;
  gf->type           = type;
  gf->score          = UNDEFDOUBLE;
  gf->range          = range;
  gf->strand         = strand;
  gf->phase          = PHASE_UNDEFINED;
  gf->attributes     = NULL;
  gf->transcripttype = TRANSCRIPT_FEATURE_TYPE_UNDETERMINED;
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
    array_add(exon_features, gf, env);
  }
  return 0;
}

void genome_feature_get_exons(GenomeFeature *gf, Array *exon_features, Env *env)
{
  int has_err;
  assert(gf && exon_features && !array_size(exon_features));
  has_err = genome_node_traverse_children((GenomeNode*) gf, exon_features,
                                          save_exon, false, env);
  assert(!has_err); /* cannot happen, because save_exon() is sane */
}

static int save_exons_and_cds(GenomeNode *gn, void *data, Env *env)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  GenomeFeature *gf;
  env_error_check(env);
  gf = (GenomeFeature*) gn;
  assert(gf && info);
  if (genome_feature_get_type(gf) == gft_exon)
    array_add(info->exon_features, gf, env);
  else if (genome_feature_get_type(gf) == gft_CDS)
    array_add(info->cds_features, gf, env);
  return 0;
}

static void set_transcript_types(Array *features)
{
  GenomeFeature *gf;
  unsigned long i;
  assert(features);
  if (array_size(features)) {
    if (array_size(features) == 1) {
      gf = *(GenomeFeature**) array_get(features, 0);
      gf->transcripttype = TRANSCRIPT_FEATURE_TYPE_SINGLE;
    }
    else {
      gf = *(GenomeFeature**) array_get(features, 0);
      gf->transcripttype = TRANSCRIPT_FEATURE_TYPE_INITIAL;
      for (i = 1; i < array_size(features) - 1; i++) {
        gf = *(GenomeFeature**) array_get(features, i);
        gf->transcripttype = TRANSCRIPT_FEATURE_TYPE_INTERNAL;
      }
      gf = *(GenomeFeature**) array_get(features, array_size(features) - 1);
      gf->transcripttype = TRANSCRIPT_FEATURE_TYPE_TERMINAL;
    }
  }
}

static int determine_transcripttypes(GenomeNode *gn, void *data, Env *env)
{
  SaveExonAndCDSInfo *info = (SaveExonAndCDSInfo*) data;
  int has_err;
  env_error_check(env);
  assert(gn && info);
  /* reset exon_features and cds_features */
  array_set_size(info->exon_features, 0);
  array_set_size(info->cds_features, 0);
  /* collect all direct children exons */
  has_err = genome_node_traverse_direct_children(gn, info, save_exons_and_cds,
                                                 env);
  assert(!has_err); /* cannot happen, because save_exon() is sane */
  /* set transcript feature type, if necessary */
  set_transcript_types(info->exon_features);
  set_transcript_types(info->cds_features);
  return 0;
}

void genome_feature_determine_transcripttypes(GenomeFeature *gf, Env *env)
{
  SaveExonAndCDSInfo info;
  int has_err;
  assert(gf);
  info.exon_features = array_new(sizeof (GenomeFeature*), env);
  info.cds_features = array_new(sizeof (GenomeFeature*), env);
  has_err = genome_node_traverse_children((GenomeNode*) gf, &info,
                                          determine_transcripttypes, false,
                                          env);
  assert(!has_err); /* cannot happen, because determine_transcripttypes() is
                       sane */
  array_delete(info.exon_features, env);
  array_delete(info.cds_features, env);
}

TranscriptFeatureType genome_feature_get_transcriptfeaturetype(GenomeFeature
                                                               *gf)
{
  assert(gf);
  return gf->transcripttype;
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
                                  const char *attr_value, Env *env)
{
  assert(gf && attr_name && attr_value);
  if (!gf->attributes)
    gf->attributes = hashtable_new(HASH_DIRECT, env_ma_free_func,
                                   env_ma_free_func, env);
  hashtable_add(gf->attributes, cstr_dup(attr_name, env),
                cstr_dup(attr_value, env), env);
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
