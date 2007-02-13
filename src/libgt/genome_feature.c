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
#include "range.h"
#include "strand.h"
#include "undef.h"
#include "xansi.h"

struct Genome_feature
{
  const GenomeNode parent_instance;
  Str *seqid,
      *source;
  GenomeFeatureType type;
  Range range;
  double score;
  Strand strand;
  Phase phase;
};

#define genome_feature_cast(GN)\
        genome_node_cast(genome_feature_class(), GN)

static void genome_feature_free(GenomeNode *gn)
{
  Genome_feature *gf = genome_feature_cast(gn);
  assert(gf);
  str_free(gf->seqid);
  str_free(gf->source);
}

static Str* genome_feature_get_seqid(GenomeNode *gn)
{
  Genome_feature *gf = genome_feature_cast(gn);
  return gf->seqid;
}

static Range genome_feature_get_range(GenomeNode *gn)
{
  Genome_feature *gf = genome_feature_cast(gn);
  return gf->range;
}

static void genome_feature_set_seqid(GenomeNode *gn, Str *seqid)
{
  Genome_feature *gf = genome_feature_cast(gn);
  assert(gf && seqid && !gf->seqid);
  gf->seqid = str_ref(seqid);
}

static void genome_feature_set_source(GenomeNode *gn, Str *source)
{
  Genome_feature *gf = genome_feature_cast(gn);
  assert(gf && source && !gf->source);
  gf->source = str_ref(source);
}

void genome_feature_set_score(Genome_feature *gf, double score)
{
  assert(gf);
  gf->score = score;
}

static void genome_feature_set_phase(GenomeNode *gn, Phase phase)
{
  Genome_feature *gf = genome_feature_cast(gn);
  assert(gf && gf->phase == PHASE_UNDEFINED);
  gf->phase = phase;
}

static void genome_feature_accept(GenomeNode *gn, GenomeVisitor *gv, Log *l)
{
  Genome_feature *gf = genome_feature_cast(gn);
  genome_visitor_visit_genome_feature(gv, gf, l);
}

const GenomeNodeClass* genome_feature_class()
{
  static const GenomeNodeClass gnc = { sizeof(Genome_feature),
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
  Genome_feature *gf = genome_feature_cast(gn);
  assert(range.start <= range.end);
  gf->seqid  = NULL;
  gf->source = NULL;
  gf->type   = type;
  gf->score  = UNDEFDOUBLE;
  gf->range  = range;
  gf->strand = strand;
  gf->phase  = PHASE_UNDEFINED;
  return gn;
}

const char* genome_feature_get_source(Genome_feature *gf)
{
  assert(gf);
  return gf->source ? str_get(gf->source) : ".";
}

GenomeFeatureType genome_feature_get_type(Genome_feature *gf)
{
  assert(gf);
  return gf->type;
}

double genome_feature_get_score(Genome_feature *gf)
{
  assert(gf);
  return gf->score;
}

Strand genome_feature_get_strand(Genome_feature *gf)
{
  assert(gf);
  return gf->strand;
}

Phase genome_feature_get_phase(Genome_feature *gf)
{
  assert(gf);
  return gf->phase;
}

static void save_exon(GenomeNode *gn, void *data)
{
  Genome_feature *gf = (Genome_feature*) gn;
  Array *exon_features = (Array*) data;
  assert(gf && exon_features);
  if (genome_feature_get_type(gf) == gft_exon) {
    array_add(exon_features, gf);
  }
}

void genome_feature_get_exons(Genome_feature *gf, Array *exon_features)
{
  assert(gf && exon_features && !array_size(exon_features));
  genome_node_traverse_children((GenomeNode*) gf, exon_features, save_exon, 0);
}

void genome_feature_set_end(Genome_feature *gf, unsigned long end)
{
  assert(gf && gf->range.start <= end);
  gf->range.end = end;
}
