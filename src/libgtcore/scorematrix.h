/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

#include "libgtcore/alpha.h"
#include "libgtcore/error.h"

typedef struct ScoreMatrix ScoreMatrix;

/* a score matrix is always defined over a given alphabet */
ScoreMatrix* scorematrix_new(Alpha*);
/* reads in a protein scorematrix from the given <path> and returns it */
ScoreMatrix* scorematrix_read_protein(const char *path, Error*);
unsigned int scorematrix_get_dimension(const ScoreMatrix*);
int          scorematrix_get_score(const ScoreMatrix*,
                                   unsigned int, unsigned int);
void         scorematrix_set_score(ScoreMatrix*,
                                   unsigned int, unsigned int, int);
void         scorematrix_show(const ScoreMatrix*, FILE*);
void         scorematrix_delete(ScoreMatrix*);

#endif
