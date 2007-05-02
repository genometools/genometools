/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <gtcore.h>
#include <libgtext/feature_visitor.h>
#include <libgtext/genome_visitor_rep.h>
#include <libgtext/sequence_region.h>

struct FeatureVisitor {
  const GenomeVisitor parent_instance;
	Hashtable *features;
};

#define feature_visitor_cast(GV)\
        genome_visitor_cast(feature_visitor_class(), GV)

static void feature_visitor_free(GenomeVisitor *gv, Env *env)
{
  FeatureVisitor *feature_visitor = feature_visitor_cast(gv);
  assert(feature_visitor);
}

/*
	Prints some information about a given genome feature object.
	Can be used as a iterator function for subtree traversal.
*/
static int print_feature_children(GenomeNode* gn, void* data, Env *env) {
	printf("%s, %lu - %lu \n", genome_feature_type_get_cstr(
																genome_feature_get_type((GenomeFeature*) gn)),
														 genome_node_get_start(gn),
														 genome_node_get_end(gn));
	/* Always returns 0 for now. */
	return 0;
}


static int feature_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                      Env *env)
{
  FeatureVisitor *v = feature_visitor_cast(gv);
	GenomeFeature* feat_from_array;
  Str *seqid;
	int has_err;
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
		array_add(hashtable_get(v->features, seqid),
						  gf, 
							env);
		
		/* Debug: Get  last added node from array, to make sure that it is the right one */
		feat_from_array = *(GenomeFeature**) array_get_last(hashtable_get(v->features, seqid));
		
		has_err = genome_node_traverse_children((GenomeNode*) feat_from_array, 
																						NULL, 
																						print_feature_children, 
																						1, 
																						env);
	}
	return has_err;
}


static int feature_visitor_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                      Env *env)
{
  FeatureVisitor *v = feature_visitor_cast(gv);
	int has_err = 0;
  env_error_check(env);
	
	if (sr == NULL) 
	{
		has_err = -1;
	}
	else
	{ 
		/* initialize new Array of subtree nodes for this sequence 
		 region and register in HT */
  	hashtable_add(v->features, 
								genome_node_get_seqid((GenomeNode*) sr), 
								array_new(sizeof(GenomeNode*),env),
								env);
	}
	return has_err;
}

const GenomeVisitorClass* feature_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (FeatureVisitor),
                                          feature_visitor_free,
                                          NULL,
                                          feature_visitor_genome_feature,
                                          feature_visitor_sequence_region,
                                          NULL };
  return &gvc;
}

GenomeVisitor* feature_visitor_new(Hashtable *features,
                                   Env *env)
{
  GenomeVisitor *gv;
  FeatureVisitor *feature_visitor;
  env_error_check(env);
  assert(features);
  gv = genome_visitor_create(feature_visitor_class(), env);
  feature_visitor = feature_visitor_cast(gv);
  feature_visitor->features = features;
  return gv;
}
