/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SCOREFUNCTION_H
#define SCOREFUNCTION_H

#include "scorematrix.h"

typedef struct ScoreFunction ScoreFunction;

ScoreFunction* scorefunction_new(ScoreMatrix*, /* takes ownership  */
                                 int deletion_score, int insertion_score);
int            scorefunction_get_score(const ScoreFunction*,
                                       unsigned char, unsigned char);
int            scorefunction_get_deletion_score(const ScoreFunction*);
int            scorefunction_get_insertion_score(const ScoreFunction*);
void           scorefunction_free(ScoreFunction*);

#endif
