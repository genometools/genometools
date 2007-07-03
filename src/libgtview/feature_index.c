/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Malte Mader <mmader@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/hashtable.h>
#include <libgtcore/range.h>
#include <libgtcore/env.h>
#include <libgtview/feature_index.h>
#include <libgtext/genome_node.h>
#include <libgtext/bsearch.h>

struct FeatureIndex
{
  Hashtable *features;
  Hashtable *ranges;
  char* firstseqid;
  unsigned int nof_sequence_regions;
};

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
  fi->features = hashtable_new(HASH_STRING, NULL, NULL, env);
  fi->ranges = hashtable_new(HASH_STRING, NULL, env_ma_free_func, env);
  return fi;
}

/*!
Prints some information about a given genome feature object.
Can be used as an iterator function for subtree traversal.
\param gn Parent GenomeNode to traverse.
\param data Pointer to memory location that could persist between calls.
\param env Pointer to Environment object.
\return Status code: 0 if ok.
*/
int genome_node_print_feature_children(GenomeNode *gn, void *data, Env *env)
{
  assert(gn);
  printf("%s, %lu - %lu \n", genome_feature_type_get_cstr(
           genome_feature_get_type((GenomeFeature*) gn)),
           genome_node_get_start(gn),
           genome_node_get_end(gn));
  /* Always returns 0 for now. */
  return 0;
}

/*!
Deletes a hashtable entry in the FeatureIndex table.
Can be used as an iterator function for hashtable_foreach().
\param key Pointer to HT key
\param value Pointer to HT value
\param data Pointer to memory location that could persist between calls.
\param env Pointer to Environment object.
\return Status code: 0 if ok.
*/
int delete_FI_row(void *key, void *value, void *data, Env *env)
{
  unsigned long i=0;
  Array *arr = (Array*) value;
  for (i = 0; i < array_size(arr); i++)
    genome_node_rec_delete(*(GenomeNode**) array_get(arr, i), env);
  array_delete(arr, env);
  /* Always returns 0 for now. */
  return 0;
}

/*!
Deletes a FeatureIndex object.
\param fi Pointer to FeatureIndex object to delete.
\param env Pointer to Environment object.
*/
void feature_index_delete(FeatureIndex *fi, Env *env)
{
  assert(fi);
  hashtable_foreach(fi->features, delete_FI_row, NULL, env);
  hashtable_delete(fi->features, env);
  hashtable_delete(fi->ranges, env);
  env_ma_free(fi, env);
}

/*!
Creates a new table entry for some sequence region.
\param fi Pointer to FeatureIndex object to add to.
\param sr Pointer to SequenceRegion object to process.
\param env Pointer to Environment object.
\return Status code: 0 if ok, -1 otherwise.
*/
int feature_index_add_sequence_region(FeatureIndex *fi,
                                      char *seqid,
                                      Env *env)
{
  assert(fi && seqid);
  int had_err = 0;
  Range *range = env_ma_malloc(env, sizeof (Range));
  env_error_check(env);
  if (fi == NULL || seqid == NULL)
  {
    had_err = -1;
  }
  else
  {
    /* initialize new Array of subtree nodes for this sequence
       region and register in HT */
    hashtable_add(fi->features,
                  seqid,
                  array_new(sizeof (GenomeNode*),env),
                  env);
    /* initialize new seqid-associated Range in a hashtable */
    hashtable_add(fi->ranges,
                  seqid,
                  range,
                  env);
    range->start = ~0UL;
    range->end = 0;
    if (fi->nof_sequence_regions++ == 0)
      fi->firstseqid = seqid;
  }
  return had_err;
}

/*!
Adds a GenomeFeature to the index, associating it with a
sequence region denoted by its identifier string.
\param fi FeatureIndex object to add to.
\param gf GenomeFeature object to add.
\param env Pointer to Environment object.
*/
void feature_index_add_genome_feature_for_seqid(FeatureIndex *fi,
                                                GenomeFeature *gf,
                                                Env *env)
{
  assert(fi && gf);
  GenomeNode *gn;
  GenomeFeature *gf_new;
  char* seqid;
  Range node_range, *seqid_range;
  env_error_check(env);

  gn = genome_node_rec_ref((GenomeNode*) gf, env);

  gf_new = (GenomeFeature*) gn;

  /* get sequence region for given GenomeFeature */
  seqid = str_get(genome_node_get_seqid(gn));

  /* create key if it is not there yet */
  if (!feature_index_has_seqid(fi, seqid, env))
    feature_index_add_sequence_region(fi, seqid, env);

  /* Add node to the appropriate array in the hashtable. */
  array_add((Array*) hashtable_get(fi->features, seqid),
            gf_new,
            env);

  /* update maximum range for given seqid */
  node_range = genome_node_get_range(gn);
  seqid_range = (Range*) hashtable_get(fi->ranges, seqid);

  if (node_range.start < seqid_range->start)
    seqid_range->start = node_range.start;
  if (node_range.end > seqid_range->end)
    seqid_range->end = node_range.end;
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
  assert(fi);
  Array* res=NULL;
  res = (Array*) hashtable_get(fi->features, seqid);
  return res;
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
  Array* base = feature_index_get_features_for_seqid(fi, seqid);
  GenomeNode* key;
  int had_err = 0, i = 0;

  assert(fi && results && seqid && (qry_range.start < qry_range.end));

  key = genome_feature_new(gft_gene, qry_range, STRAND_UNKNOWN, NULL,
                           UNDEF_ULONG, env);

 for (i = 0; i < array_size(base); i++) {
   GenomeNode *gn = *(GenomeNode**) array_get(base, i);
   Range r = genome_node_get_range(gn);
   if (range_overlap(r, qry_range)) {
     array_add(results, gn, env);
   }
 }

 genome_node_delete(key, env);

 return had_err;
}

/*
Outputs a row in a hashtable, with the key being a
cstring and the value an Array.
Output will be "<key> -> <number of elements in array>".
Meant to be used in hashtable_foreach().
*/
static int print_index_row(void *key, void *value, void *data, Env *env)
{
  int i;
  printf("%s -> %lu\n", (char*) key, array_size((Array*) value));
  for (i=0;i<array_size((Array*) value);i++)
  {
    GenomeNode* gn = *(GenomeNode**) array_get((Array*) value, i);
    Range r = genome_node_get_range(gn);
    printf("%s, %s: %lu - %lu \n",
           (char*) key,
            genome_feature_type_get_cstr(
                   genome_feature_get_type((GenomeFeature* )gn)),
            r.start, r.end);
  }
  return 0;
}

/*!
Simple output of array counts for all sequence regions in a FI.
\param fi FeatureIndex object to lookup in.
*/
void feature_index_print_contents(FeatureIndex *fi, Env *env)
{
  assert(fi);
  hashtable_foreach(fi->features, print_index_row, NULL, env);
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
  assert(fi && seqid);
  return (hashtable_get(fi->features, seqid) != NULL);
}

/*!
Returns the first sequence region identifier added to the index.
\param fi FeatureIndex object to lookup in.
\return sequence region identifier cstr pointer
*/
char* feature_index_get_first_seqid(FeatureIndex *fi)
{
  assert(fi);
  return (fi->firstseqid ? fi->firstseqid : "");
}

Range feature_index_get_range_for_seqid(FeatureIndex *fi,
                                        char *seqid)
{
  assert(fi && seqid);
  Range ret, *range;
  ret.start = 1;
  ret.end = ~0UL;

  range = (Range*) hashtable_get(fi->ranges, seqid);
  if (range)
  {
    ret.start = range->start;
    ret.end = range->end;
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
  Range r1, r2, r3, r4, r5, check_range;
  Str *seqid1, *seqid2;
  int had_err=0;

  /* Generating some ranges */
  r1.start=100; r1.end=1000;
  r2.start=100; r2.end=300;
  r3.start=500; r3.end=1000;
  r4.start=600; r4.end=1200;
  r5.start=600; r5.end=1000;

  /* Generating sequnce ids as c-strings */
  seqid1 = str_new_cstr("test1", env);
  seqid2 = str_new_cstr("test2", env);

  /* Generating a new genome_feature with the property gft_gene and the range r1 ... */
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

  /*Add a sequence region the fetaure index and test if it has really been added*/
  feature_index_add_sequence_region(fi, "test1", env);
  ensure(had_err, feature_index_has_seqid(fi, "test1", env));

  feature_index_add_sequence_region(fi, "test2", env);
  ensure(had_err, feature_index_has_seqid(fi, "test2", env));

  /*Tests if we get a empty data structure for every added sequence region*/
  ensure(had_err, feature_index_get_features_for_seqid(fi, "test1"));
  ensure(had_err, feature_index_get_features_for_seqid(fi, "test2"));
  ensure(had_err, array_size(feature_index_get_features_for_seqid(fi, "test1")) == 0);
  ensure(had_err, array_size(feature_index_get_features_for_seqid(fi, "test2")) == 0);

  /*Add features to every sequence region and test if the according datastructures are not empty anymore. As we have added one genome_feature to every sequence region the size has to be one. */
  feature_index_add_genome_feature_for_seqid(fi, (GenomeFeature*) gn1, env);
  ensure(had_err, array_size(feature_index_get_features_for_seqid(fi, "test1")) == 1);

  feature_index_add_genome_feature_for_seqid(fi, (GenomeFeature*) gn2, env);

  ensure(had_err, array_size(feature_index_get_features_for_seqid(fi, "test2")) == 1);

  /* test feature_index_get_first_seqid() */
  ensure(had_err, feature_index_get_first_seqid(fi));
  ensure(had_err, strcmp("test1", feature_index_get_first_seqid(fi)) == 0);

  check_range = feature_index_get_range_for_seqid(fi, "test1");
  ensure(had_err, check_range.start == 100 && check_range.end == 1000);

  /*delete all generated objects*/
  feature_index_delete(fi, env);
  genome_node_rec_delete(gn1, env);
  genome_node_rec_delete(gn2, env);
  str_delete(seqid1, env);
  str_delete(seqid2, env);
  return had_err;
}
