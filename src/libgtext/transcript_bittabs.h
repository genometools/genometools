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

#ifndef TRANSCRIPT_BITTABS_H
#define TRANSCRIPT_BITTABS_H

#include "libgtcore/bittab.h"

/* a container class for transcript bittabs */
typedef struct TranscriptBittabs TranscriptBittabs;

/* create an empy container */
TranscriptBittabs* transcript_bittabs_new(unsigned long size_all,
                                          unsigned long size_single,
                                          unsigned long size_initial,
                                          unsigned long size_internal,
                                          unsigned long size_terminal);

/* return the bittab for all exons */
Bittab*            transcript_bittabs_get_all(const TranscriptBittabs*);

/* return the bittab for single exons */
Bittab*            transcript_bittabs_get_single(const TranscriptBittabs*);

/* return the bittab for initial exons */
Bittab*            transcript_bittabs_get_initial(const TranscriptBittabs*);

/* return the bittab for internal exons */
Bittab*            transcript_bittabs_get_internal(const TranscriptBittabs*);

/* return the bittab for terminal exons */
Bittab*            transcript_bittabs_get_terminal(const TranscriptBittabs*);

void               transcript_bittabs_delete(TranscriptBittabs*);

#endif
