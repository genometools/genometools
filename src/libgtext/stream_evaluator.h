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

#ifndef STREAM_EVALUATOR_H
#define STREAM_EVALUATOR_H

#include "libgtext/genome_stream.h"

typedef struct StreamEvaluator StreamEvaluator;

StreamEvaluator* stream_evaluator_new(GenomeStream *reality,
                                      GenomeStream *prediction, bool nuceval,
                                      bool evalLTR, unsigned long LTRdelta);
/* if <gv> is not NULL, it visits all nodes from reality and the prediction */
int              stream_evaluator_evaluate(StreamEvaluator*, bool verbose,
                                           bool exondiff, GenomeVisitor *gv,
                                           Error*);
void             stream_evaluator_show(StreamEvaluator*, FILE*);
void             stream_evaluator_delete(StreamEvaluator*);

#endif
