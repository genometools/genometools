/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
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
			  GenomeNode* parent,
			  Env* env);
Str* track_get_title(Track* track);
Line* get_next_free_line(Track* track,
                         Range r,
			 Env* env);
Array* track_get_lines(Track* track);
int track_get_number_of_lines(Track* track);
void track_finish(Track* track,
                  Env* env);
void track_delete(Track* track,
                  Env* env);
void print_track(Track* track);
int track_unit_test(Env* env);

#endif

