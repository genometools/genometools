/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef STREAM_EVALUATOR_H
#define STREAM_EVALUATOR_H

#include "libgtext/genome_stream.h"

typedef struct StreamEvaluator StreamEvaluator;

StreamEvaluator* stream_evaluator_new(GenomeStream *reality,
                                      GenomeStream *prediction,
                                      bool evalLTR, unsigned long LTRdelta,
                                      Env*);
/* if <gv> is not NULL, it visits all nodes from reality and the prediction */
int              stream_evaluator_evaluate(StreamEvaluator*, bool verbose,
                                           bool exondiff, GenomeVisitor *gv,
                                           Env*);
void             stream_evaluator_show(StreamEvaluator*, FILE*);
void             stream_evaluator_delete(StreamEvaluator*, Env*);

#endif
