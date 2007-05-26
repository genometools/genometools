/*
   Copyright (c) 2007 Sascha Steinbiss, Christin Schaerfer, Malte Mader
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

  line = get_next_free_line(track->lines, gn);

  if(line == NULL)
  {
    line = line_new(env);
    array_add(track->lines, line, env);
  }
  
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
Line* get_next_free_line(Array *lines, GenomeNode *gn)
{
  int i;

  for(i=0; i<array_size(lines); i++)
  {
    if(!line_is_occupied((Line*) array_get(lines, i), gn))
    {
      return (Line*) array_get(lines, i);
    }
  }
  return NULL;
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
    line_delete((Line*) array_get(track->lines, i));
  }

  array_delete(track->lines, env);
  env_ma_free(track, env);
}

