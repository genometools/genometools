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


/*!
Creates a new FeatureIndex object. 
\param env Pointer to Environment object.
\return Pointer to a new FeatureIndex object.
*/
FeatureIndex* feature_index_new(Env *env) 
{
  FeatureIndex *fi;
  env_error_check(env);
  fi = env_ma_malloc(env, sizeof(FeatureIndex));
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
static int print_feature_children(GenomeNode* gn, void* data, Env *env) 
{
  assert(gn);
  env_error_check(env);
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
  env_error_check(env);
  for(i=0;i<array_size(arr);i++)
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
  env_error_check(env);
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
int feature_index_add_sequence_region(FeatureIndex *fi, SequenceRegion *sr, Env *env)
{
  assert(fi && sr);
  int has_err = 0;
  env_error_check(env);
  if (fi == NULL || sr == NULL) 
  {
    has_err = -1;
  }
  else
  { 
    /* initialize new Array of subtree nodes for this sequence 
       region and register in HT */
    hashtable_add(fi->features, 
         str_get(genome_node_get_seqid((GenomeNode*) sr)), 
         array_new(sizeof(GenomeNode*),env),
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
void feature_index_add_genome_feature_for_seqid(FeatureIndex *fi, GenomeFeature* gf, Env *env)
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
  range1 = genome_node_get_range((GenomeNode*) gn1);
  range2 = genome_node_get_range(*(GenomeNode**) gn2);
  if(range_overlap(range1 ,range2))
    return 0;
  if(range1.end < range2.start)
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
int feature_index_get_features_for_range(FeatureIndex *fi, Array* results, char* seqid, Range qry_range, Env* env)
{
 Array* base = feature_index_get_features_for_seqid(fi, seqid);
 GenomeNode* key;
 int has_err = 0;
 
assert(fi && results && seqid && (qry_range.start < qry_range.end)); 
 
 key = genome_feature_new(gft_gene, qry_range, STRAND_UNKNOWN,
                                          NULL, UNDEF_ULONG, env);
                                          
 bsearch_all(results, 
             key, 
             array_get_space(base), 
             array_size(base), 
             sizeof(GenomeNode*),
             compare_for_overlap, 
             env);
 
 genome_node_delete(key, env);
 return has_err;
}



/* 
Outputs a row in a hashtable, with the key being a cstring and the value an Array.
Output will be "<key> -> <number of elements in array>".
Meant to be used in hashtable_foreach().
*/
static int print_index_row(void* key, void* value, void* data, Env* env)
{
  printf("%s -> %lu\n", (char*) key, array_size((Array*) value));
  return 0;
}

/*
Simple output of hashtable counts for debugging reasons.
*/
void feature_index_print_contents(FeatureIndex *fi, Env* env)
{
  assert(fi);
  hashtable_foreach(fi->features, print_index_row, NULL, env);
}
