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

#ifndef TRANSCRIPT_EXONS_H
#define TRANSCRIPT_EXONS_H

#include "extended/transcript_bittabs.h"
#include "extended/transcript_counts.h"

/* a container class for transcript exon arrays */
typedef struct GtTranscriptExons GtTranscriptExons;

GtTranscriptExons* gt_transcript_exons_new(void);

/* return the exon array for all exons */
GtArray*           gt_transcript_exons_get_all(const GtTranscriptExons*);

/* return the exon array for single exons */
GtArray*           gt_transcript_exons_get_single(const GtTranscriptExons*);

/* return the exon array for initial exons */
GtArray*           gt_transcript_exons_get_initial(const GtTranscriptExons*);

/* return the exon array for internal exons */
GtArray*           gt_transcript_exons_get_internal(const GtTranscriptExons*);

/* return the exon array for terminal exons */
GtArray*           gt_transcript_exons_get_terminal(const GtTranscriptExons*);

void               gt_transcript_exons_sort(const GtTranscriptExons*);

GtTranscriptCounts* gt_transcript_exons_uniq_in_place_count(GtTranscriptExons*);

bool               gt_transcript_exons_are_sorted(const GtTranscriptExons*);
GtTranscriptBittabs*
                   gt_transcript_exons_create_bittabs(const GtTranscriptExons*);

void               gt_transcript_exons_delete(GtTranscriptExons*);

#endif
