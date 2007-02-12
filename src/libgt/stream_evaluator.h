/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef STREAM_EVALUATOR_H
#define STREAM_EVALUATOR_H

#include "genome_stream.h"

typedef struct Stream_evaluator Stream_evaluator;

/* takes ownership of the given streams */
Stream_evaluator* stream_evaluator_new(Genome_stream *reality,
                                       Genome_stream *prediction);
void              stream_evaluator_evaluate(Stream_evaluator*,
                                            bool verbose, bool exondiff);
void              stream_evaluator_show(Stream_evaluator*, FILE*);
void              stream_evaluator_free(Stream_evaluator*);

#endif
