/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TRANSCRIPT_EVALUATORS_H
#define TRANSCRIPT_EVALUATORS_H

#include <libgtext/evaluator.h>
#include <libgtext/transcript_exons.h>

/* a container class for transcript evaluators */
typedef struct TranscriptEvaluators TranscriptEvaluators;

TranscriptEvaluators* transcript_evaluators_new(Env*);

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

void                  transcript_evaluators_delete(TranscriptEvaluators*, Env*);

#endif
