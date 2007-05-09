#include <libgtcore/hashtable.h>
#include <libgtcore/range.h>
#include <libgtext/feature_index.h>
#include <libgtext/genome_node.h>
#include <libgtext/bsearch.h>

struct FeatureIndex
{
  Hashtable *features;
};

/*
Creates a new FeatureIndex object. 
*/
FeatureIndex* feature_index_new(Env *env) 
{
	FeatureIndex *fi;
	env_error_check(env);
	fi = env_ma_malloc(env, sizeof(FeatureIndex));
	fi->features = hashtable_new(HASH_STRING, NULL, NULL, env);
	return fi;
}


/*
Prints some information about a given genome feature object.
Can be used as an iterator function for subtree traversal.
*/
static int print_feature_children(GenomeNode* gn, void* data, Env *env) {
	printf("%s, %lu - %lu \n", genome_feature_type_get_cstr(
					 genome_feature_get_type((GenomeFeature*) gn)),
				   genome_node_get_start(gn),
				   genome_node_get_end(gn));
	/* Always returns 0 for now. */
	return 0;
}

int delete_rows(void* key, void* value, void* data, Env* env)
{
  unsigned long i=0;
  Array *arr = (Array*) value;
  for(i=0;i<array_size(arr);i++)
  {
  	genome_node_rec_delete(*(GenomeNode**) array_get(arr, i), env);
  }
  array_delete(arr, env);
  return 0;
}

/*
Deletes a FeatureIndex object.
*/
void feature_index_delete(FeatureIndex *fi, Env *env) 
{
  hashtable_foreach(fi->features, delete_rows, NULL, env);
	hashtable_delete(fi->features, env);
	env_ma_free(fi, env);
}

/*
Creates a new empty array in the feature index for a new SequenceRegion.
Returns 0 if ok, otherwise <0.
*/
int feature_index_add_sequence_region(FeatureIndex *fi, SequenceRegion *sr, Env *env)
{
	int has_err = 0;
	if (sr == NULL) 
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

/*
Adds a new GenomeFeature object to the feature index. The sequence region is determined from the 
entry itself. 
*/
void feature_index_add_genome_feature_for_seqid(FeatureIndex *fi, GenomeFeature* gf, Env *env)
{
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

/*
Returns an array of GenomeNodes for a given sequence region identifier.
*/
Array* feature_index_get_features_for_seqid(FeatureIndex* fi, char* seqid)
{
	return (Array*) hashtable_get(fi->features, seqid);	
}


static int compare_for_overlap(const void* gn1, const void* gn2)
{
  Range range1, range2;
  range1 = genome_node_get_range((GenomeNode*) gn1);
  range2 = genome_node_get_range(*(GenomeNode**) gn2);
  if(range_overlap(range1 ,range2))
    return 0;
  if(range1.end < range2.start)
    return -1;
  return 1;
}


/*

*/
int feature_index_get_features_for_range(FeatureIndex *fi, Array* results, char* seqid, unsigned long range_start, unsigned long range_end, Env* env)
{
 Array* base = feature_index_get_features_for_seqid(fi, seqid);
 GenomeNode* key;
 int has_err = 0;
 
 Range qry_range, tmp_range;
 qry_range.start = range_start;
 qry_range.end = range_end;
 
 key = genome_feature_new(gft_gene, qry_range, STRAND_UNKNOWN,
                                          NULL, UNDEF_ULONG, env);
 tmp_range = genome_node_get_range(key);                                        
 printf("tmp start: %lu, end: %lu\n", tmp_range.start, tmp_range.end);
 
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
Simple output of hashtable contents for debugging reasons.
*/
void feature_index_print_contents(FeatureIndex *fi, Env* env)
{
  hashtable_foreach(fi->features, print_index_row, NULL, env);
}
