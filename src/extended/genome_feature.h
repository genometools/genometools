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

#ifndef GENOME_FEATURE_H
#define GENOME_FEATURE_H

/* implements the ``genome node'' interface */
typedef struct GT_GenomeFeature GT_GenomeFeature;

#include "core/range.h"
#include "core/phase.h"
#include "core/strand.h"
#include "core/strarray.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"
#include "extended/transcript_feature_type.h"

typedef void (*AttributeIterFunc)(const char *attr_name, const char *attr_value,
                                  void *data);

const GT_GenomeNodeClass* gt_genome_feature_class(void);
GT_GenomeNode*        gt_genome_feature_new(GT_Str *seqid, const char *type,
                                            GT_Range, GT_Strand);
GT_GenomeNode*        gt_genome_feature_new_pseudo(GT_GenomeFeature*);
/* Return the ``standard gene'' (mainly for testing purposes). */
GT_GenomeNode*        gt_genome_feature_new_standard_gene(void);
const char*           gt_genome_feature_get_source(GT_GenomeFeature*);
const char*           gt_genome_feature_get_attribute(GT_GenomeNode *gn,
                                                      const char *attr_name);
/* Return a GT_StrArray containing the used attribute names. */
GT_StrArray*          gt_genome_feature_get_attribute_list(GT_GenomeFeature*);
const char*           gt_genome_feature_get_type(GT_GenomeFeature*);
bool                  gt_genome_feature_has_type(GT_GenomeFeature*,
                                                 const char*);
bool                  gt_genome_feature_score_is_defined(const
                                                         GT_GenomeFeature*);
bool                  gt_genome_feature_is_multi(const GT_GenomeFeature*);
bool                  gt_genome_feature_is_pseudo(const GT_GenomeFeature*);
void                  gt_genome_feature_make_multi_representative(const
                                                             GT_GenomeFeature*);
void                  gt_genome_feature_set_multi_representative(
                                                             GT_GenomeFeature*,
                                                             GT_GenomeFeature*);
GT_GenomeFeature*     gt_genome_feature_get_multi_representative(
                                                             GT_GenomeFeature*);
float                 gt_genome_feature_get_score(GT_GenomeFeature*);
GT_Strand             gt_genome_feature_get_strand(GT_GenomeFeature*);
Phase                 gt_genome_feature_get_phase(GT_GenomeFeature*);
void                  gt_genome_feature_get_exons(GT_GenomeFeature*,
                                                  GT_Array *exon_features);
void                  gt_genome_feature_determine_transcripttypes(
                                                             GT_GenomeFeature*);
TranscriptFeatureType gt_genome_feature_get_transcriptfeaturetype(
                                                             GT_GenomeFeature*);
void                  gt_genome_feature_set_source(GT_GenomeNode*,
                                                   GT_Str *source);
void                  gt_genome_feature_set_phase(GT_GenomeNode*, Phase);
void                  gt_genome_feature_set_end(GT_GenomeFeature*,
                                                unsigned long);
void                  gt_genome_feature_set_score(GT_GenomeFeature*, float);
void                  gt_genome_feature_unset_score(GT_GenomeFeature*);
void                  gt_genome_feature_add_attribute(GT_GenomeFeature*,
                                                      const char *attr_name,
                                                      const char *attr_value);
void                  gt_genome_feature_foreach_attribute(GT_GenomeFeature*,
                                                          AttributeIterFunc,
                                                          void *data);
bool                  gt_genome_feature_has_CDS(const GT_GenomeFeature*);
bool                  gt_genome_feature_has_splice_site(const
                                                        GT_GenomeFeature*);
double                gt_genome_feature_average_splice_site_prob(const
                                                             GT_GenomeFeature*);
/* Returns true, if the given features have the same seqid, feature type, range,
   strand, and phase. */
bool                  gt_genome_features_are_similar(GT_GenomeFeature*,
                                                     GT_GenomeFeature*);
int                   gt_genome_feature_unit_test(GT_Error*);

#endif
