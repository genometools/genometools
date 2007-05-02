#include <libgtcore/hashtable.h>
#include <libgtext/feature_index.h>
#include <libgtext/genome_node.h>

struct FeatureIndex
{
  Hashtable *features;
};

FeatureIndex* feature_index_new(Env *env) 
{
	FeatureIndex *fi;
	env_error_check(env);
	fi = env_ma_malloc(env, sizeof(FeatureIndex));
	fi->features = hashtable_new(HASH_STRING, NULL, (FreeFunc) array_delete, env);
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

void feature_index_delete(FeatureIndex *fi, Env *env) 
{
	hashtable_delete(fi->features, env);
	env_ma_free(fi, env);
}

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
								genome_node_get_seqid((GenomeNode*) sr), 
								array_new(sizeof(GenomeNode*),env),
								env);
	}
	
	return has_err;
}

int feature_index_add_genome_feature_for_seqid(FeatureIndex *fi, Str* seqid, GenomeFeature* gf, Env *env)
{
	GenomeFeature* feat_from_array;
	int has_err = 0;
	env_error_check(env);

	/* get sequence region for given GenomeFeature */
	seqid = genome_node_get_seqid((GenomeNode*) gf);
	if(seqid == NULL) 
	{
		has_err = -1; 
	} 
	else 
	{
		/* Add node to the appropriate array in the hashtable. */
		array_add(hashtable_get(fi->features, seqid),
						  gf, 
							env);
		
		/* Debug: Get  last added node from array, to make sure that it is the right one */
		feat_from_array = *(GenomeFeature**) array_get_last(hashtable_get(fi->features, seqid));
		/* Debug: print children of current node */
		genome_node_traverse_children((GenomeNode*) feat_from_array, 
																						NULL, 
																						print_feature_children, 
																						1, 
																						env);
	}
	return has_err;
}

int feature_index_get_features_for_seqid(FeatureIndex* fi, Array* feature_array, Str* seqid)
{
	int has_err = 0;
	if(seqid == NULL || fi == NULL) 
	{
		has_err = -1; 
	} 
	else 
	{
		feature_array = hashtable_get(fi->features, seqid);
	}	
	return  has_err;
}
