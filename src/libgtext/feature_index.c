/*
  Copyright (c) 2007 Sascha Steinbiss, Christin Schaerfer, Malte Mader
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/hashtable.h>
#include <libgtcore/range.h>
#include <libgtext/feature_index.h>
#include <libgtext/genome_node.h>
#include <libgtext/bsearch.h>

struct FeatureIndex
{
  Hashtable *features;
};

typedef struct
{
  char* searchkey;
  bool found;
} ht_search_params;

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
  fi->features = hashtable_new(HASH_STRING, NULL, NULL, env);
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
int genome_node_print_feature_children(GenomeNode* gn, void* data, Env *env)
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
int delete_FI_row(void* key, void* value, void* data, Env* env)
{
  unsigned long i=0;
  Array *arr = (Array*) value;
  for (i=0;i<array_size(arr);i++)
  {
  	genome_node_rec_delete(*(GenomeNode**) array_get(arr, i), env);
  }
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
                                      char* seqid,
                                      Env *env)
{
  assert(fi && seqid);
  int has_err = 0;
  env_error_check(env);
  if (fi == NULL || seqid == NULL)
  {
    has_err = -1;
  }
  else
  {
    /* initialize new Array of subtree nodes for this sequence
       region and register in HT */
    hashtable_add(fi->features,
         seqid,
         array_new(sizeof (GenomeNode*),env),
         env);
  }
  return has_err;
}

/*!
Adds a GenomeFeature to the index, associating it with a
sequence region denoted by its identifier string.
\param fi FeatureIndex object to add to.
\param gf GenomeFeature object to add.
\param env Pointer to Environment object.
*/
void feature_index_add_genome_feature_for_seqid(FeatureIndex *fi,
                                                GenomeFeature* gf,
                                                Env *env)
{
  assert(fi && gf);
  GenomeNode *gn;
  GenomeFeature *gf_new;
  char* seqid;
  env_error_check(env);

  gn = genome_node_rec_ref((GenomeNode*) gf, env);

  gf_new = (GenomeFeature*) gn;

  /* get sequence region for given GenomeFeature */
  seqid = str_get(genome_node_get_seqid(gn));

  /* create key if it is not there yet */
  if (!feature_index_has_seqid(fi, seqid, env))
	{
	  feature_index_add_sequence_region(fi, seqid, env);
	}

  /* Add node to the appropriate array in the hashtable. */
  array_add(hashtable_get(fi->features, seqid),
            gf_new,
            env);
}

/*!
Returns an array of GenomeFeatures associated with a given
sequence region identifier.
\param fi FeatureIndex object to add to.
\param seqid Sequence region identifier to lookup.
\return Pointer to the result array.
*/
Array* feature_index_get_features_for_seqid(FeatureIndex* fi, char* seqid)
{
  assert(fi);
  Array* res=NULL;
  res = (Array*) hashtable_get(fi->features, seqid);
  return res;
}

/*!
Comparator for GenomeNodes. Overlaps are treated as equality.
\param gn1 Pointer to GenomeNode object.
\param gn1 Pointer to GenomeNode object.
\return 0 if nodes overlap, -1 if gn1
ends strictly left of gn2, 1 if gn2 starts strictly right of gn1.
*/
static int compare_for_overlap(const void* gn1, const void* gn2)
{
  assert(gn1 && gn2);
  Range range1, range2;
  range1 = genome_node_get_range(*(GenomeNode**) gn1);
  range2 = genome_node_get_range(*(GenomeNode**) gn2);
  if (range_overlap(range1 ,range2))
    return 0;
  if (range1.end < range2.start)
    return -1;
  return 1;
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
                                         Array* results,
                                         char* seqid,
                                         Range qry_range,
                                         Env* env)
{
 Array* base = feature_index_get_features_for_seqid(fi, seqid);
 GenomeNode* key;
 int has_err = 0;

assert(fi && results && seqid && (qry_range.start < qry_range.end));

 key = genome_feature_new(gft_gene, qry_range, STRAND_UNKNOWN,
                                          NULL, UNDEF_ULONG, env);

 bsearch_all(results,
             &key,
             array_get_space(base),
             array_size(base),
             sizeof (GenomeFeature*),
             compare_for_overlap,
             env);

 genome_node_delete(key, env);
 return has_err;
}

/*
Outputs a row in a hashtable, with the key being a
cstring and the value an Array.
Output will be "<key> -> <number of elements in array>".
Meant to be used in hashtable_foreach().
*/
static int print_index_row(void* key, void* value, void* data, Env* env)
{
  printf("%s -> %lu\n", (char*) key, array_size((Array*) value));
  return 0;
}

/*!
Simple output of array counts for all sequence regions in a FI.
\param fi FeatureIndex object to lookup in.
*/
void feature_index_print_contents(FeatureIndex *fi, Env* env)
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
bool feature_index_has_seqid(FeatureIndex* fi, char* seqid, Env* env)
{
  assert(fi && seqid);

  return (hashtable_get(fi->features, seqid) != NULL);
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
	Range r1, r2, r3, r4, r5;
	Str *seqid1, *seqid2;
	int has_err=0;
	
	/* Generating some ranges */
	r1.start=100; r1.end=1000;
	r2.start=100; r2.end=300;
	r3.start=500; r3.end=1000;
	r4.start=600; r4.end=1200;
	r5.start=600; r5.end=1000;
	
	/* Generating sequnce ids as c-strings */
	seqid1 = str_new_cstr("test1", env);
	seqid2 = str_new_cstr("test2", env);
		
	/* Generating a new genome_feature with the property gft_gene an the range r1 ... */
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

  ensure(has_err, fi);
  ensure(has_err, !feature_index_has_seqid(fi, "test1", env));
  ensure(has_err, !feature_index_has_seqid(fi, "test2", env));
	
	/*Add a sequence region the fetaure index and test if it has really been added*/
	feature_index_add_sequence_region(fi, "test1", env);
  ensure(has_err, feature_index_has_seqid(fi, "test1", env));

	feature_index_add_sequence_region(fi, "test2", env);
  ensure(has_err, feature_index_has_seqid(fi, "test2", env));
  
  /*Tests if we get a empty data structure for every added sequence region*/
	ensure(has_err, feature_index_get_features_for_seqid(fi, "test1"));
	ensure(has_err, feature_index_get_features_for_seqid(fi, "test2"));
	ensure(has_err, array_size(feature_index_get_features_for_seqid(fi,
                                                                  "test1")
																																	) == 0);
	ensure(has_err, array_size(feature_index_get_features_for_seqid(fi,
	                                                                "test2")
																																	) == 0);
	
	/*Add features to every sequence region and test if the according 
	datastructures are not empty anymore. As we have added one 
	genome_feature to every sequence region the size has to be one.
	*/
	feature_index_add_genome_feature_for_seqid(fi, (GenomeFeature*) gn1, env);
	ensure(has_err, array_size(feature_index_get_features_for_seqid(fi,
	                                                                "test1")
																																	) == 1);

	feature_index_add_genome_feature_for_seqid(fi, (GenomeFeature*) gn2, env);

	ensure(has_err, array_size(feature_index_get_features_for_seqid(fi,
	                                                                "test2")
																																	) == 1);
	
	/*delete all generated objects*/
	feature_index_delete(fi, env);
	genome_node_rec_delete(gn1, env);
	genome_node_rec_delete(gn2, env);
	str_delete(seqid1, env);
	str_delete(seqid2, env);
  return has_err;
}
