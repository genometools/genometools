/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TRANSCRIPT_COUNTS_H
#define TRANSCRIPT_COUNTS_H

#include <gtcore.h>

/* a container class for transcript count arrays */
typedef struct TranscriptCounts TranscriptCounts;

/* create an empy container */
TranscriptCounts* transcript_counts_new(Env*);

/* return the count array for all exons */
Array*            transcript_counts_get_all(const TranscriptCounts*);
void              transcript_counts_set_all(TranscriptCounts*, Array*);

/* return the count array for single exons */
Array*            transcript_counts_get_single(const TranscriptCounts*);
void              transcript_counts_set_single(TranscriptCounts*, Array*);

/* return the count array for initial exons */
Array*            transcript_counts_get_initial(const TranscriptCounts*);
void              transcript_counts_set_initial(TranscriptCounts*, Array*);

/* return the count array for internal exons */
Array*            transcript_counts_get_internal(const TranscriptCounts*);
void              transcript_counts_set_internal(TranscriptCounts*, Array*);

/* return the count array for terminal exons */
Array*            transcript_counts_get_terminal(const TranscriptCounts*);
void              transcript_counts_set_terminal(TranscriptCounts*, Array*);

void              transcript_counts_delete(TranscriptCounts*, Env*);

#endif
