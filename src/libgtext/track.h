/*
   Copyright (c) Sascha Steinbiss, Malte Mader, Christin Schaerfer
   Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TRACK_H
#define TRACK_H

#include <libgtcore/str.h>
#include <libgtcore/array.h>
#include <libgtext/genome_node.h>
#include <libgtext/config.h>
#include <libgtext/line.h>

typedef struct Track Track;

Track* track_new(Str* title,
                 Env* env);
void track_insert_element(Track* track, 
                          GenomeNode* gn, 
			  Config* cfg, 
			  Env* env);
Str* track_get_title(Track* track);
Line* get_next_free_line(Track* track,
                         GenomeNode *gn);
Array* track_get_lines(Track* track);
void track_delete(Track* track,
                  Env* env);

#endif

