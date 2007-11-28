/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef TRANSCRIPT_EVALUATORS_H
#define TRANSCRIPT_EVALUATORS_H

#include "libgtext/evaluator.h"
#include "libgtext/transcript_exons.h"

/* a container class for transcript evaluators */
typedef struct TranscriptEvaluators TranscriptEvaluators;

TranscriptEvaluators* transcript_evaluators_new(void);

/* return the evaluator for all exons */
Evaluator*            transcript_evaluators_get_all(const
                                                    TranscriptEvaluators*);

/* return the evaluator for single exons */
Evaluator*            transcript_evaluators_get_single(const
                                                       TranscriptEvaluators*);

/* return the evaluator for initial exons */
Evaluator*            transcript_evaluators_get_initial(const
                                                        TranscriptEvaluators*);

/* return the evaluator for internal exons */
Evaluator*            transcript_evaluators_get_internal(const
                                                         TranscriptEvaluators*);

/* return the evaluator for terminal exons */
Evaluator*            transcript_evaluators_get_terminal(const
                                                         TranscriptEvaluators*);

void                  transcript_evaluators_add_actuals(const
                                                        TranscriptEvaluators*,
                                                        const TranscriptExons*);

void                  transcript_evaluators_delete(TranscriptEvaluators*);

#endif
