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

#include "libgtext/transcript_bittabs.h"
#include "libgtext/transcript_counts.h"

/* a container class for transcript exon arrays */
typedef struct TranscriptExons TranscriptExons;

TranscriptExons*   transcript_exons_new(void);

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

TranscriptCounts*  transcript_exons_uniq_in_place_count(TranscriptExons*);

bool               transcript_exons_are_sorted(const TranscriptExons*);
TranscriptBittabs* transcript_exons_create_bittabs(const TranscriptExons*);

void               transcript_exons_delete(TranscriptExons*);

#endif
