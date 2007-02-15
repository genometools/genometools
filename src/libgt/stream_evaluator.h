/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef STREAM_EVALUATOR_H
#define STREAM_EVALUATOR_H

#include "genome_stream.h"

typedef struct StreamEvaluator StreamEvaluator;

/* takes ownership of the given streams */
StreamEvaluator* stream_evaluator_new(GenomeStream *reality,
                                      GenomeStream *prediction);
int              stream_evaluator_evaluate(StreamEvaluator*, bool verbose,
                                           bool exondiff, Error*);
void             stream_evaluator_show(StreamEvaluator*, FILE*);
void             stream_evaluator_free(StreamEvaluator*);

#endif
