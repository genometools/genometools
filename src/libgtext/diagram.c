/*
  Copyright (c) 2007 Sascha Steinbiss, Christin Schaerfer, Malte Mader
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/diagram.h>
#include <libgtext/track.h>
#include <libgtcore/hashtable.h>
#include <libgtext/genome_node.h>
#include <libgtcore/str.h>
#include <libgtext/genome_feature_type.h>
#include <libgtext/genome_feature.h>

struct Diagram
{
  Hashtable *tracks;
  Config *config;
  Range range;
};

/*!
Check if a genome node is in the given Range of the diagram. If yes, hand it
over to the corresponding track object. If the track doesn't exist, generate
it.
\param pointer to the diagram object.
\param pointer to the array of genome nodes.
\param env Pointer to Environment object.
*/
static int visit_child(GenomeNode* gn, void* diagram, Env* env)
{
  const char* feature_type;
  Diagram* d = (Diagram*) diagram;
  GenomeFeatureType type;
  GenomeFeature* gf = (GenomeFeature*) gn;
  Track* track;
  Range gn_range;
  Str* track_type;
  
  /*fetch the type of the given genome node*/
  type = genome_feature_get_type(gf);
  feature_type = genome_feature_type_get_cstr(type);
  
  gn_range = genome_node_get_range(gn);
  
  /*Check the Range of the genome node. If the start of the gn is before the 
  start of the diagram or the end of the gn is after the end of the 
  digram, set the corresponding value of the gn node to the value of the
  diagram.*/
  if(gn_range.end < d->range.start || gn_range.start > d->range.end)
  {
    return 0;
  }
  else if (gn_range.start < d->range.start)
  {
    gn_range.start = d->range.start;
  }
  else if (gn_range.end > d->range.end)
  {
    gn_range.end = d->range.end;
  }
  
  /*deliver the genome node to the track with the corresponding type.
  If the track doesn't exit, generate it.*/
  if (hashtable_get(d->tracks, feature_type) == NULL)
  {
    track_type = str_new_cstr((char*) feature_type, env);
    track = track_new(track_type, env);
    track_insert_element(track, gn, d->config, env);
  }
  else
  {	
    track = hashtable_get(d->tracks, feature_type);
    track_insert_element(track, gn, d->config, env);
  }
  
  return 0;
}

/*!
Iterating through the array of genome nodes and traversing the genome node
trees. Calling the function genome node with the current genome node.
\param pointer to the diagram object.
\param pointer to the array of genome nodes.
\param env Pointer to Environment object.
*/
static void diagram_build(Diagram* diagram, Array* features, Env* env)
{
  int i=0;
  for(i=0;i<array_size(features);i++)
  {
    genome_node_traverse_children(**(GenomeNode***) array_get(features, i),                                  
				  diagram,
                                  visit_child,
				  true,
				  env);
  }
}


/*
Initialize a new diagram object.
\param pointer to the array of genome nodes.
\param the given range of the diagam.
\param pointer to the configuration object.
\param env Pointer to Environment object.
*/
Diagram* diagram_new(Array* features, Range range, Config* config, Env* env)
{
	
  Diagram *diagram;
  env_error_check(env);
  diagram = env_ma_malloc(env, sizeof (Diagram));
  diagram->tracks = hashtable_new(HASH_STRING, NULL, NULL, env);
  diagram->range = range;
  diagram_build(diagram, features, env);
  return diagram; 
}

/*
Update the configuration object with new settings.
\param pointer to the diagram object.
\param pointer to the configuration object.
\param env Pointer to Environment object.
*/
void diagram_set_config(Diagram* diagram, Config* config, Env* env)
{
  diagram->config = config;
}

/*
Delete the diagramm object.
\param pointer to the diagram object.
\param env Pointer to Environment object.
*/
void diagram_delete(Diagram* diagram, Env* env)
{
  hashtable_delete(diagram->tracks, env);
  env_ma_free(diagram, env);
}
