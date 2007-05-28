/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/track.h>

struct Track
{
  Str *title;
  Array *lines;
};

/*!
Creates a new Track object.
\param title Title of the Track object.
\param env Pointer to Environment object.
\return Pointer to a new Track object.
*/
Track* track_new(Str *title,
                 Env *env)
{
  Track *track;
  env_error_check(env);
  track = env_ma_malloc(env, sizeof (Track));
  track->title = title;
  track->lines = array_new(sizeof (Line*), env);
  return track;
}

/*!
Inserts an element into a Track object
\param track Track in which to insert element 
\param gn Pointer to GenomeNode object
\param cfg Pointer to Config file
\param env Pointer to Environment object
*/
void track_insert_element(Track *track,
                          GenomeNode *gn,
			  Config *cfg,
			  Env *env)
{
  Line *line;

  line = get_next_free_line(track, gn, env);
  line_insert_element(line, gn, cfg, env);
}

/*!
Returns Track title
\param track Pointer to Track object
\teturn Pointer to Title String object
*/
Str* track_get_title(Track *track)
{
  return track->title;
}

/*!
Gets the next unoccupied Line object
\param lines Array with Line objects
\param gn Pointer to GenomeNode object
\return Pointer to unoccupied Line object
*/
Line* get_next_free_line(Track* track, GenomeNode *gn, Env* env)
{
  int i;
  Line* line;

  for(i=0; i<array_size(track->lines); i++)
  {
    line = *(Line**) array_get(track->lines, i);
    if(!line_is_occupied(line, gn))
    {
      return line;
    }
  }
  line = line_new(env);
  array_add(track->lines, line, env);
  return line;
}

/*!
Returns Array with Pointer to Line objects
\param track Pointer to Track object
\return Pointer to Array
*/
Array* track_get_lines(Track* track)
{
  return track->lines;
}


/*!
Delets Track
\param Track Pointer to Track object to delet
\param env Pointer to Environment object
*/
void track_delete(Track *track,
                  Env *env)
{
  int i;

  if(!track) return;

  for(i=0; i<array_size(track->lines); i++)
  {
    line_delete(*(Line**) array_get(track->lines, i), env);
  }

  array_delete(track->lines, env);
  env_ma_free(track, env);
}

/*!
Prints all Lines of a Track object
uparam track Pointer to Track object to print
*/
void print_track(Track* track)
{ 
  int i;
    
  for(i=0; i<array_size(track->lines); i++)
  { 
    print_line(*(Line**) array_get(track->lines, i));
    printf("\n");
   }
}

/*!
 * Tests Track Class
 * \param env Pointer to Environment object
 * */
int track_unit_test(Env* env)
{
  Range r1, r2, r3; 
  Array* lines;
  int has_err = 0;

  r1.start = 10;
  r1.end = 50;

  r2.start = 60;
  r2.end = 80;

  r3.start = 70;
  r3.end = 100;

  GenomeNode* gn1 = genome_feature_new(gft_exon, r1, STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn2 = genome_feature_new(gft_intron, r2, STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn3 = genome_feature_new(gft_intron, r3, STRAND_FORWARD, NULL, 0, env);
			                    
  Line* l1 = line_new(env);
  Line* l2 = line_new(env);

  Str* title = str_new(env);
  str_set(title, "exon", env);
  Str* s = str_new(env);
  str_set(s, "foo", env);

  Track* t = track_new(title, env);
  Track* t2 = track_new(s, env);

  /* test track_insert_elements 
     (implicit test of get_next_free_line) */
  ensure(has_err, (0 == array_size(track_get_lines(t))));
  track_insert_element(t, gn1, NULL, env);
  ensure(has_err, (1 == array_size(track_get_lines(t))));
  track_insert_element(t, gn2, NULL, env);
  ensure(has_err, (1 == array_size(track_get_lines(t))));
  track_insert_element(t, gn3, NULL, env);
  ensure(has_err, (2 == array_size(track_get_lines(t))));

  /* test track_get_title */
  ensure(has_err, (0 == str_cmp(title, track_get_title(t))));
  ensure(has_err, !(0 == str_cmp(s, track_get_title(t))));

  /* test track_get_lines */
  lines = track_get_lines(t2);
  ensure(has_err, (0 == array_size(lines)));
  track_insert_element(t2, gn1, NULL, env);
  lines = track_get_lines(t2);
  ensure(has_err, (1 == array_size(lines)));
  lines = track_get_lines(t);
  ensure(has_err, (2 == array_size(lines)));

  line_delete(l1, env);
  line_delete(l2, env);
  track_delete(t, env);
  track_delete(t2, env);
  str_delete(title, env);
  str_delete(s, env);
  genome_node_delete(gn1, env);
  genome_node_delete(gn2, env);
  genome_node_delete(gn3, env);

  return has_err;
}

