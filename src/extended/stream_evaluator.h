/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "extended/node_stream_api.h"
#include "extended/node_visitor_api.h"

typedef struct GtStreamEvaluator GtStreamEvaluator;

GtStreamEvaluator* gt_stream_evaluator_new(GtNodeStream *reference,
                                           GtNodeStream *prediction,
                                           bool nuceval, bool evalLTR,
                                           GtUword LTRdelta);
/* if <nv> is not NULL, it visits all nodes from reference and the prediction */
int                gt_stream_evaluator_evaluate(GtStreamEvaluator*,
                                                bool verbose, bool exondiff,
                                                bool exondiffcollapsed,
                                                GtNodeVisitor *nv, GtError*);
void               gt_stream_evaluator_show(GtStreamEvaluator*, GtFile*);
void               gt_stream_evaluator_delete(GtStreamEvaluator*);

#endif
