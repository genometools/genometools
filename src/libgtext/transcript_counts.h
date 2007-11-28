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

#ifndef TRANSCRIPT_COUNTS_H
#define TRANSCRIPT_COUNTS_H

#include "libgtcore/array.h"

/* a container class for transcript count arrays */
typedef struct TranscriptCounts TranscriptCounts;

/* create an empy container */
TranscriptCounts* transcript_counts_new(void);

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

void              transcript_counts_delete(TranscriptCounts*);

#endif
