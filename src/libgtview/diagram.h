/*
  Copyright (c) 2007 Malte Mader <mmader@zbh.uni-hamburg.de>
  Copyright (c) 2007 Sascha Steinbiss, Christin Schaerfer, Malte Mader
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef DIAGRAM_H
#define DIAGRAM_H

#include "libgtview/config.h"
#include "libgtview/block.h"
#include "libgtview/feature_index.h"
#include "libgtcore/array.h"
#include "libgtcore/range.h"
#include "libgtcore/hashtable.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_feature_type.h"

typedef struct Diagram Diagram;

/* Create a new diagram object representing the genome nodes in <features>
   overlapping with <range>. */
Diagram*    diagram_new(Array *features, Range range, Config*, Env*);
Range       diagram_get_range(Diagram*);
void        diagram_set_config(Diagram*, Config*, Env*);
Hashtable*  diagram_get_tracks(const Diagram*);
int         diagram_get_total_lines(const Diagram*, Env*);
int         diagram_get_number_of_tracks(const Diagram*);
int         diagram_unit_test(Env*);
void        diagram_delete(Diagram*, Env*);

#endif
