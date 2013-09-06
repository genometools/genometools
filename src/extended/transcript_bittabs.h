/*
  Copyright (c) 2007 Gordon Gremme <gordon@gremme.org>
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

#ifndef TRANSCRIPT_BITTABS_H
#define TRANSCRIPT_BITTABS_H

#include "core/bittab.h"

/* a container class for transcript bittabs */
typedef struct GtTranscriptBittabs GtTranscriptBittabs;

/* create an empy container */
GtTranscriptBittabs* gt_transcript_bittabs_new(GtUword size_all,
                                               GtUword size_single,
                                               GtUword size_initial,
                                               GtUword size_internal,
                                               GtUword size_terminal);

/* return the bittab for all exons */
GtBittab*            gt_transcript_bittabs_get_all(const GtTranscriptBittabs*);

/* return the bittab for single exons */
GtBittab*            gt_transcript_bittabs_get_single(const
                                                      GtTranscriptBittabs*);

/* return the bittab for initial exons */
GtBittab*            gt_transcript_bittabs_get_initial(const
                                                       GtTranscriptBittabs*);

/* return the bittab for internal exons */
GtBittab*            gt_transcript_bittabs_get_internal(const
                                                        GtTranscriptBittabs*);

/* return the bittab for terminal exons */
GtBittab*            gt_transcript_bittabs_get_terminal(const
                                                        GtTranscriptBittabs*);

void                 gt_transcript_bittabs_delete(GtTranscriptBittabs*);

#endif
