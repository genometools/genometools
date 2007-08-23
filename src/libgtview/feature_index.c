/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>,
  Copyright (c) 2007 Malte Mader <mmader@stud.zbh.uni-hamburg.de>,
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/
/**
 * \if INTERNAL \file feature_index.c \endif
 * \author Malte Mader <mmader@zbh.uni-hamburg.de>
 * \author Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
 * \author Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
 */

#include <string.h>
#include "libgtcore/ensure.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/minmax.h"
#include "libgtcore/range.h"
#include "libgtcore/undef.h"
#include "libgtview/feature_index.h"
#include "libgtext/genome_node.h"
#include "libgtext/bsearch.h"

struct FeatureIndex
{
  Hashtable *regions;
  char* firstseqid;
  unsigned int nof_sequence_regions;
};

typedef struct
{
  Array *features;
  SequenceRegion *region;
  Range dyn_range;
} RegionInfo;

/*!
Deletes a hashtable entry in the FeatureIndex table.
Can be used as a FreeFunc.
*/
int region_info_delete(void *data, Env *env)
{
  unsigned long i=0;
  RegionInfo *info = (RegionInfo*) data;
  for (i = 0; i < array_size(info->features); i++)
    genome_node_rec_delete(*(GenomeNode**) array_get(info->features, i), env);
  array_delete(info->features, env);
  genome_node_rec_delete((GenomeNode*)info->region, env);
  env_ma_free_func(info, env);
  return 0;
}

/*!
Creates a new FeatureIndex object.
\param env Pointer to Environment object.
\return Pointer to a new FeatureIndex object.
*/
FeatureIndex* feature_index_new(Env *env)
{
  FeatureIndex *fi;
  env_error_check(env);
  fi = env_ma_malloc(env, sizeof (FeatureIndex));
  fi->nof_sequence_regions = 0;
  fi->firstseqid = NULL;
  fi->regions  = hashtable_new(HASH_STRING, NULL,
                               (FreeFunc) region_info_delete, env);
  return fi;
}

/*!
Deletes a FeatureIndex object.
\param fi Pointer to FeatureIndex object to delete.
\param env Pointer to Environment object.
*/
void feature_index_delete(FeatureIndex *fi, Env *env)
{
  assert(fi != NULL);
  hashtable_delete(fi->regions, env);
  env_ma_free(fi, env);
}

/*!
Creates a new table entry for some sequence region.
\param fi Pointer to FeatureIndex object to add to.
\param sr Pointer to SequenceRegion object to process.
\param env Pointer to Environment object.
*/
void feature_index_add_sequence_region(FeatureIndex *fi,
                                      SequenceRegion *sr,
                                      Env *env)
{
  char *seqid;
  RegionInfo *info;

  assert(fi != NULL && sr != NULL);

  seqid = str_get(genome_node_get_seqid((GenomeNode*) sr));
  info = env_ma_malloc(env, sizeof (RegionInfo));
  info->region = (SequenceRegion*) genome_node_rec_ref((GenomeNode*) sr, env);
  info->features = array_new(sizeof (GenomeNode*),env);
  info->dyn_range.start = ~0UL;
  info->dyn_range.end   = 0;
  hashtable_add(fi->regions,
                seqid,
                info,
                env);
  if (fi->nof_sequence_regions++ == 0)
    fi->firstseqid = seqid;
}

/*!
Adds a GenomeFeature to the index, associating it with a
sequence region denoted by its identifier string.
\param fi FeatureIndex object to add to.
\param gf GenomeFeature object to add.
\param env Pointer to Environment object.
*/
void feature_index_add_genome_feature(FeatureIndex *fi,
                                      GenomeFeature *gf,
                                      Env *env)
{
  GenomeNode *gn;
  GenomeFeature *gf_new;
  char* seqid;
  Range node_range;
  RegionInfo *info;
  env_error_check(env);

  assert(fi != NULL && gf != NULL);

  gn = genome_node_rec_ref((GenomeNode*) gf, env);
  gf_new = (GenomeFeature*) gn;
  /* get information about seqid and range */
  node_range = genome_node_get_range(gn);
  seqid = str_get(genome_node_get_seqid(gn));
  /* entry for the seqid must already exist */
  assert(feature_index_has_seqid(fi, seqid, env));
  info = (RegionInfo*) hashtable_get(fi->regions, seqid);
  /* add node to the appropriate array in the hashtable. */
  array_add(info->features,
            gf_new,
            env);
  /* update dynamic range */
  info->dyn_range.start = MIN(info->dyn_range.start, node_range.start);
  info->dyn_range.end = MAX(info->dyn_range.end, node_range.end);
}

/*!
Returns an array of GenomeFeatures associated with a given
sequence region identifier.
\param fi FeatureIndex object to add to.
\param seqid Sequence region identifier to lookup.
\return Pointer to the result array.
*/
Array* feature_index_get_features_for_seqid(FeatureIndex *fi, char *seqid)
{
  RegionInfo *res;

  assert(fi != NULL);

  res = (RegionInfo*) hashtable_get(fi->regions, seqid);
  return (res != NULL ? res->features : NULL);
}

/*!
Looks up relevant GenomeFeatures in a given range inside of a given
sequence region.
\param fi FeatureIndex object to lookup in.
\param results Array object to store results in.
\param seqid Sequence region identifier to extract features for.
\param qry_range Query range.
\return 0 if ok.
*/
int feature_index_get_features_for_range(FeatureIndex *fi,
                                         Array *results,
                                         char *seqid,
                                         Range qry_range,
                                         Env *env)
{
  Array* base;
  GenomeNode* key;
  unsigned long i = 0;

  assert(fi != NULL && results != NULL);

  base = feature_index_get_features_for_seqid(fi, seqid);
  if (!base) return -1;
  assert(fi && results && seqid && base && (qry_range.start < qry_range.end));
  key = genome_feature_new(gft_gene, qry_range, STRAND_UNKNOWN, NULL,
                           UNDEF_ULONG, env);
  for (i = 0; i < array_size(base); i++)
  {
    GenomeNode *gn = *(GenomeNode**) array_get(base, i);
    Range r = genome_node_get_range(gn);
    if (range_overlap(r, qry_range))
      array_add(results, gn, env);
  }
  genome_node_delete(key, env);
  return 0;
}

/*!
Checks whether a given cstr is a key in the FeatureIndex.
\param fi FeatureIndex object to lookup in.
\param seqid Sequence region identifier to check.
\param env Pointer to Environment object.
\return TRUE if found, FALSE otherwise.
*/
bool feature_index_has_seqid(FeatureIndex *fi, char *seqid, Env *env)
{
  assert(fi != NULL);
  return (hashtable_get(fi->regions, seqid) != NULL);
}

/*!
Returns the first sequence region identifier added to the index.
\param fi FeatureIndex object to lookup in.
\return sequence region identifier cstr pointer
*/
char* feature_index_get_first_seqid(FeatureIndex *fi)
{
  assert(fi != NULL);
  return fi->firstseqid;
}

Range feature_index_get_range_for_seqid(FeatureIndex *fi,
                                        char *seqid)
{
  Range ret;
  RegionInfo *info;

  assert(fi != NULL);

  info = (RegionInfo*) hashtable_get(fi->regions, seqid);
  if (info && (info->dyn_range.start != ~0UL
            && info->dyn_range.end != 0))
  {
    ret.start = info->dyn_range.start;
    ret.end = info->dyn_range.end;
  }
  else
  {
    return genome_node_get_range((GenomeNode*) info->region);
  }
  return ret;
}

/*!
tests all public functions from feature_index.
\param env Pointer to Environment object.
\return 0 if ok.
*/
int feature_index_unit_test(Env* env)
{
  /* First we have to create some objects that we can use for testing. */
  GenomeNode *gn1, *gn2, *ex1, *ex2, *ex3, *cds1;
  FeatureIndex *fi;
  Range r1, r2, r3, r4, r5, check_range, rs;
  Str *seqid1, *seqid2;
  SequenceRegion *sr1, *sr2;
  int had_err=0;

  /* Generating some ranges */
  r1.start=100UL; r1.end=1000UL;
  r2.start=100UL; r2.end=300UL;
  r3.start=500UL; r3.end=1000UL;
  r4.start=600UL; r4.end=1200UL;
  r5.start=600UL; r5.end=1000UL;
  rs.start=100UL; rs.end=1200UL;

  /* Generating sequnce ids as c-strings */
  seqid1 = str_new_cstr("test1", env);
  seqid2 = str_new_cstr("test2", env);

  sr1 = (SequenceRegion*) sequence_region_new(seqid1, rs, NULL, 0, env);
  sr2 = (SequenceRegion*) sequence_region_new(seqid2, rs, NULL, 0, env);

  /* Generating a new genome_feature with the property gft_gene and the range r1
     ... */
  gn1 = genome_feature_new(gft_gene, r1, STRAND_UNKNOWN,
                                            NULL, UNDEF_ULONG, env);
  /* ... and assign a sequence id to the new genome_feature-object. */
  genome_node_set_seqid((GenomeNode*) gn1, seqid1);

  gn2 = genome_feature_new(gft_gene, r4, STRAND_UNKNOWN,
                                            NULL, UNDEF_ULONG, env);
  genome_node_set_seqid((GenomeNode*) gn2, seqid2);

  ex1 = genome_feature_new(gft_exon, r2, STRAND_UNKNOWN,
                                            NULL, UNDEF_ULONG, env);
  genome_node_set_seqid((GenomeNode*) ex1, seqid1);

  ex2 = genome_feature_new(gft_exon, r3, STRAND_UNKNOWN,
                                            NULL, UNDEF_ULONG, env);
  genome_node_set_seqid((GenomeNode*) ex2, seqid1);

  ex3 = genome_feature_new(gft_exon, r4, STRAND_UNKNOWN,
                                            NULL, UNDEF_ULONG, env);
  genome_node_set_seqid((GenomeNode*) ex3, seqid2);

  cds1 = genome_feature_new(gft_CDS, r5, STRAND_UNKNOWN,
                                            NULL, UNDEF_ULONG, env);
  genome_node_set_seqid((GenomeNode*) cds1, seqid2);

  /* Determine the structure of our feature tree */
  genome_node_is_part_of_genome_node(gn1, ex1, env);
  genome_node_is_part_of_genome_node(gn1, ex2, env);
  genome_node_is_part_of_genome_node(gn2, ex3, env);
  genome_node_is_part_of_genome_node(gn2, cds1, env);

  /*Create a new feature index on which we can perfom some tests*/
  fi  = feature_index_new(env);

  ensure(had_err, fi);
  ensure(had_err, !feature_index_has_seqid(fi, "test1", env));
  ensure(had_err, !feature_index_has_seqid(fi, "test2", env));

  /* Add a sequence region to the feature index and test if it has really been
     added */
  feature_index_add_sequence_region(fi, sr1, env);
  ensure(had_err, feature_index_has_seqid(fi, "test1", env));

  feature_index_add_sequence_region(fi, sr2, env);
  ensure(had_err, feature_index_has_seqid(fi, "test2", env));

  /* Tests if we get a empty data structure for every added sequence region*/
  ensure(had_err, feature_index_get_features_for_seqid(fi, "test1"));
  ensure(had_err, feature_index_get_features_for_seqid(fi, "test2"));
  ensure(had_err,
         array_size(feature_index_get_features_for_seqid(fi, "test1")) == 0);
  ensure(had_err,
         array_size(feature_index_get_features_for_seqid(fi, "test2")) == 0);

  /* Add features to every sequence region and test if the according
     datastructures are not empty anymore. As we have added one genome_feature
     to every sequence region the size has to be one. */
  feature_index_add_genome_feature(fi, (GenomeFeature*) gn1, env);
  ensure(had_err,
         array_size(feature_index_get_features_for_seqid(fi, "test1")) == 1UL);

  feature_index_add_genome_feature(fi, (GenomeFeature*) gn2, env);

  ensure(had_err,
         array_size(feature_index_get_features_for_seqid(fi, "test2")) == 1UL);

  /* test feature_index_get_first_seqid() */
  ensure(had_err, feature_index_get_first_seqid(fi));
  ensure(had_err, strcmp("test1", feature_index_get_first_seqid(fi)) == 0);

  check_range = feature_index_get_range_for_seqid(fi, "test1");
  ensure(had_err, check_range.start == 100UL && check_range.end == 1000UL);

  ensure(had_err, feature_index_get_features_for_seqid(fi, "test1") != NULL);
  ensure(had_err, feature_index_get_features_for_seqid(fi, "noexist") == NULL);

  /* delete all generated objects */
  feature_index_delete(fi, env);
  genome_node_rec_delete(gn1, env);
  genome_node_rec_delete(gn2, env);
  genome_node_rec_delete((GenomeNode*) sr1, env);
  genome_node_rec_delete((GenomeNode*) sr2, env);
  str_delete(seqid1, env);
  str_delete(seqid2, env);
  return had_err;
}
