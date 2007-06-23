/*
  Copyright (c) 2007 Malte Mader <mmader@zbh.uni-hamburg.de>
                     Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
                     Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/diagram.h>
#include <libgtext/track.h>
#include <libgtext/genome_node.h>
#include <libgtcore/str.h>
#include <libgtext/genome_feature_type.h>
#include <libgtext/genome_feature.h>
#include <libgtext/feature_index.h>

struct Diagram
{
  Hashtable *tracks;
  int nof_tracks;
  Config *config;
  Range range;
};

typedef struct
  {
    GenomeNode *parent;
    Diagram *diagram;

  } GenomeNodeChildren;

Range diagram_get_range(Diagram* diagram)
{
  assert(diagram);
  return diagram->range;
}

/*
Fetching the number of lines from a track object and add number
 of lines to data.
*/
static int diagram_add_tracklines(void *key, void *value, void *data,
                                  Env* env)
{
  int* add;
  add = data;
  *add += track_get_number_of_lines(value);
  return 0;
}

/*
deletes every given track object.
*/
static int diagram_track_delete(void *key, void *value, void *data, Env* env)
{
  track_delete((Track*) value, env);
  return 0;
}

/*!
Insert a genome node into a track. If the track to a
 corresponding genome type does not exist, create a new track.
\param gn pointer to a genome node object.
\param parent pointer the parent of the genome node object.
\param d pointer to the diagram object.
\param env Pointer to Environment object.
*/
static void insert_genome_node_into_track(GenomeNode* gn,
                                          GenomeNode* parent,
                                          Diagram* d,
                                          Env* env)
{
  const char* feature_type;
  GenomeFeatureType type;
  GenomeFeature* gf = (GenomeFeature*) gn;
  Track* track;
  Str* track_type;
  Range elem_range;

  /* discard elements that do not overlap with visible range */
  elem_range = genome_node_get_range(gn);
  if (!range_overlap(d->range, elem_range)) return;

  /* fetch the type of the given genome node */
  type = genome_feature_get_type(gf);
  feature_type = genome_feature_type_get_cstr(type);

  /* coming soon:
  if ( hide_always(feature_type) ||
      (type_range(feature_type) < d->range.start - d->range.start) ) return;
  */

  if (config_cstr_in_list(d->config,"collapse","to_parent",feature_type,
                             env))
  {
    type = genome_feature_get_type((GenomeFeature*) parent);
    feature_type = genome_feature_type_get_cstr(type);
  }

  /* deliver the genome node to the track with the corresponding type.
   If the track does not exist, generate it.
   */
  if (hashtable_get(d->tracks, feature_type) == NULL)
  {
    track_type = str_new_cstr((char*) feature_type, env);
    track = track_new(track_type, env);
    track_insert_element(track, gn, d->config, parent, env);
    hashtable_add(d->tracks, (char*) feature_type, track, env);
    d->nof_tracks++;
    str_delete(track_type,env);
  }
  else
  {
    track = hashtable_get(d->tracks, feature_type);
    track_insert_element(track, gn, d->config, parent, env);
  }
}

/*!
Insert a genome node into a track. If the genome node has children,
insert them with a recursiv call of visit_child.
\param gn pointer to the diagram object.
\param genome_node_children struct, containing pointers
       to a genome node and a diagram.
\param env Pointer to Environment object.
*/
static int visit_child(GenomeNode* gn, void* genome_node_children, Env* env)
{

  GenomeNodeChildren* genome_node_info;
  genome_node_info = (GenomeNodeChildren*) genome_node_children;
  
  if (genome_node_has_children(gn))
  {
    GenomeNode *oldparent = genome_node_info->parent;
    insert_genome_node_into_track(gn,
                                  genome_node_info->parent,
                                  genome_node_info->diagram,
                                  env);
    genome_node_info->parent = gn;
    genome_node_traverse_direct_children(gn,
                                         genome_node_info,
                                         visit_child,
                                         env);
    genome_node_info->parent = oldparent;
  }
  else
  {
  insert_genome_node_into_track(gn,
                                genome_node_info->parent,
                                genome_node_info->diagram,
                                env);
  }
  return 0;

}

/*
creates Lines for each track, iterator function
*/
static int finish_track(void *key, void *value, void *data, Env* env)
{
  track_finish((Track*) value, env);
  return 0;
}

/*!
Travering a genome node tree with depth first search.
 Furthermore the parent of a genome node is known.
\param gn pointer to the diagram object.
\param genome_node_children struct, containing pointers
 to a genome node and a diagram.
\param env Pointer to Environment object.
*/
static int traverse_genome_nodes(GenomeNode* gn,
                                 void* genome_node_children,
				 Env* env)
{
  assert(genome_node_children);
  GenomeNodeChildren* genome_node_info;
  genome_node_info = (GenomeNodeChildren*) genome_node_children;
  genome_node_info->parent = gn;

  if (genome_node_has_children(gn))
  {
    genome_node_traverse_direct_children(gn,
                                         genome_node_info,
                                         visit_child,
                                         env);
  }
  else
  {
    return 0;
  }
  return 0;
}

/*!
Iterating through the array of genome nodes and traversing the genome node
trees and insert always the first root node into a track.
\param pointer to the diagram object.
\param pointer to the array of genome nodes.
\param env Pointer to Environment object.
*/
static void diagram_build(Diagram* diagram, Array* features, Env* env)
{
  GenomeNodeChildren genome_node_children;
  genome_node_children.diagram = diagram;
  int i=0;
  for (i=0;i<array_size(features);i++)
  {

    GenomeNode *current_root = *(GenomeNode**) array_get(features,i);
    insert_genome_node_into_track(current_root,
                                  NULL,
                                  diagram,
                                  env);

    traverse_genome_nodes(current_root, &genome_node_children, env);
  }
	/* create lines for each track in diagram */
  hashtable_foreach(diagram->tracks, finish_track, NULL, env);
}

/*
Initialize a new diagram object.
\param features pointer to the array of genome nodes.
\param range the given range of the diagam.
\param config pointer to the configuration object.
\param env Pointer to Environment object.
*/
Diagram* diagram_new(Array* features, Range range, Config* config, Env* env)
{
  Diagram *diagram;
  env_error_check(env);
  diagram = env_ma_malloc(env, sizeof (Diagram));
  diagram->tracks = hashtable_new(HASH_STRING, NULL, NULL, env);
  diagram->nof_tracks = 0;
  diagram->config = config;
  diagram->range = range;
  diagram_build(diagram, features, env);
  return diagram;
}

/*
Update the configuration object with new settings.
\param diagram pointer to the diagram object.
\param config pointer to the configuration object.
\param env Pointer to Environment object.
*/
void diagram_set_config(Diagram* diagram, Config* config, Env* env)
{
  assert(diagram && config && env);
  diagram->config = config;
}

/*
Delivers the hashtable with the stored tracks.
\param diagram pointer to the diagram object.
*/
Hashtable* diagram_get_tracks(Diagram* diagram)
{
  assert(diagram);
  return diagram->tracks;
}

/*
Returns the number of all lines in the diagram.
\param diagram pointer to the diagram object.
\param env Pointer to Environment object.
*/
int diagram_get_total_lines(Diagram* diagram, Env* env)
{
  assert(diagram && env);
  int total_lines;
  total_lines = 0;
  hashtable_foreach(diagram->tracks, diagram_add_tracklines, &total_lines,
                    env);
  return total_lines;
}

/*
Returns the number of all lines in the diagram.
\param diagram pointer to the diagram object.
\param env Pointer to Environment object.
*/
int diagram_get_number_of_tracks(Diagram* diagram)
{
  assert(diagram);
  return diagram->nof_tracks;
}

/*
Delete the diagram object.
\param diagram pointer to the diagram object.
\param env Pointer to Environment object.
*/
void diagram_delete(Diagram* diagram, Env* env)
{
  assert(diagram && env);
  hashtable_foreach(diagram->tracks, diagram_track_delete, NULL, env);
  hashtable_delete(diagram->tracks, env);
  env_ma_free(diagram, env);
}

/*
generate a feature index test structure and test the diagram functions.
*/
int diagram_unit_test(Env* env)
{

  GenomeNode *gn1, *gn2, *ex1, *ex2, *ex3, *cds1;
  FeatureIndex *fi;
  Range r1, r2, r3, r4, r5, dr1;
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

  /* Generating a new genome_feature with the property gft_gene and
   the range r1 ... */
  gn1 = genome_feature_new(gft_gene, r1, STRAND_UNKNOWN, NULL, UNDEF_ULONG,
                           env);

  /* ... and assign a sequence id to the new genome_feature-object. */
  genome_node_set_seqid((GenomeNode*) gn1, seqid1);

  gn2 = genome_feature_new(gft_gene, r4, STRAND_UNKNOWN, NULL, UNDEF_ULONG,
                           env);

  genome_node_set_seqid((GenomeNode*) gn2, seqid2);

  ex1 = genome_feature_new(gft_exon, r2, STRAND_UNKNOWN, NULL, UNDEF_ULONG,
                           env);

  genome_node_set_seqid((GenomeNode*) ex1, seqid1);

  ex2 = genome_feature_new(gft_exon, r3, STRAND_UNKNOWN, NULL, UNDEF_ULONG,
                           env);

  genome_node_set_seqid((GenomeNode*) ex2, seqid1);

  ex3 = genome_feature_new(gft_exon, r4, STRAND_UNKNOWN, NULL, UNDEF_ULONG,
                           env);

  genome_node_set_seqid((GenomeNode*) ex3, seqid2);

  cds1 = genome_feature_new(gft_CDS, r5, STRAND_UNKNOWN, NULL, UNDEF_ULONG,
                            env);

  genome_node_set_seqid((GenomeNode*) cds1, seqid2);

  /* Determine the structure of our feature tree */
  genome_node_is_part_of_genome_node(gn1, ex1, env);
  genome_node_is_part_of_genome_node(gn1, ex2, env);
  genome_node_is_part_of_genome_node(gn2, ex3, env);
  genome_node_is_part_of_genome_node(gn2, cds1, env);

  /*Create a new feature index on which we can perfom some tests*/
  fi  = feature_index_new(env);

  /*Add a sequence region the feature index and test if
   it has really been added*/
  feature_index_add_sequence_region(fi, "test1", env);

  feature_index_add_sequence_region(fi, "test2", env);

  /*Add features to every sequence region.*/
  feature_index_add_genome_feature_for_seqid(fi, (GenomeFeature*) gn1, env);
  feature_index_add_genome_feature_for_seqid(fi, (GenomeFeature*) gn2, env);

  /* Generating the Range for the diagram*/
  dr1.start=400; dr1.end=900;

  /*get the features for the test1 sequence region*/
  Array *features;
  features = array_new(sizeof (GenomeFeature*), env);
  feature_index_get_features_for_range(fi,features,"test1",dr1,env);

  /*create a config object*/
  Config *cfg;
  Str *luafile = str_new_cstr("config.lua",env);
  cfg = config_new(env, false);
  config_load_file(cfg, luafile, env);

  /*create a diagram object and test it*/
  Diagram *dia;
  dia = diagram_new(features,dr1,cfg,env);

  ensure(has_err, dia->config != NULL);
  ensure(has_err, dia->range.start == 400);
  ensure(has_err, dia->range.end == 900);
  if (!config_cstr_in_list(dia->config,"collapse","to_parent","gene",
                             env))
  {
  ensure(has_err, hashtable_get(dia->tracks,"gene") != NULL);
  }

  if (!config_cstr_in_list(dia->config,"collapse","to_parent","exon",
                             env))
  {
  ensure(has_err, hashtable_get(dia->tracks,"exon") != NULL);
  }
  ensure(has_err, range_compare(diagram_get_range(dia),dr1) == 0);

  /*get the features for the test2 sequence region*/
  Array *features2;
  features2 = array_new(sizeof (GenomeFeature*), env);
  feature_index_get_features_for_range(fi,features2,"test2",dr1,env);

  /*create a diagram object and test it*/
  Diagram *dia2;
  dia2 = diagram_new(features2,dr1,cfg,env);
  ensure(has_err, dia->range.start == 400);
  ensure(has_err, dia->range.end == 900);

  if (!config_cstr_in_list(dia2->config,"collapse","to_parent","gene",
                             env))
  {
  ensure(has_err, hashtable_get(dia2->tracks,"gene") != NULL);
  }

  if (!config_cstr_in_list(dia2->config,"collapse","to_parent","exon",
                             env))
  {
  ensure(has_err, hashtable_get(dia2->tracks,"exon") != NULL);
  }

  if (!config_cstr_in_list(dia2->config,"collapse","to_parent","CDS",
                             env))

  {
  ensure(has_err, hashtable_get(dia2->tracks,"CDS") != NULL);
  }
  ensure(has_err, range_compare(diagram_get_range(dia),dr1) == 0)

  /*delete all generated objects*/
  str_delete(luafile, env);
  config_delete(cfg, env);
  diagram_delete(dia,env);
  diagram_delete(dia2,env);
  feature_index_delete(fi, env);
  array_delete(features,env);
  array_delete(features2,env);
  genome_node_rec_delete(gn1, env);
  genome_node_rec_delete(gn2, env);
  str_delete(seqid1, env);
  str_delete(seqid2, env);
return has_err;
}
