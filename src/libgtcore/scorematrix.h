/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

#include <libgtcore/alpha.h>
#include <libgtcore/env.h>

typedef struct ScoreMatrix ScoreMatrix;

/* a score matrix is always defined over a given alphabet */
ScoreMatrix* scorematrix_new(Alpha*, Env*);
/* reads in a protein scorematrix from the given path and returns it */
ScoreMatrix* scorematrix_read_protein(const char *path, Env*);
int          scorematrix_get_score(const ScoreMatrix*,
                                   unsigned char, unsigned char);
void         scorematrix_set_score(ScoreMatrix*,
                                   unsigned char, unsigned char, int);
void         scorematrix_show(const ScoreMatrix*, FILE*);
void         scorematrix_delete(ScoreMatrix*, Env*);

#endif
