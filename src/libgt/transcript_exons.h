/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TRANSCRIPT_EXONS_H
#define TRANSCRIPT_EXONS_H

#include "array.h"
#include "transcript_bittabs.h"
#include "transcript_counts.h"

/* a container class for transcript exon arrays */
typedef struct TranscriptExons TranscriptExons;

TranscriptExons*   transcript_exons_new(Env*);

/* return the exon array for all exons */
Array*             transcript_exons_get_all(const TranscriptExons*);

/* return the exon array for single exons */
Array*             transcript_exons_get_single(const TranscriptExons*);

/* return the exon array for initial exons */
Array*             transcript_exons_get_initial(const TranscriptExons*);

/* return the exon array for internal exons */
Array*             transcript_exons_get_internal(const TranscriptExons*);

/* return the exon array for terminal exons */
Array*             transcript_exons_get_terminal(const TranscriptExons*);

void               transcript_exons_sort(const TranscriptExons*);

TranscriptCounts*  transcript_exons_uniq_in_place_count(TranscriptExons*, Env*);

bool               transcript_exons_are_sorted(const TranscriptExons*);
TranscriptBittabs* transcript_exons_create_bittabs(const TranscriptExons*,
                                                   Env*);

void               transcript_exons_delete(TranscriptExons*, Env*);

#endif
