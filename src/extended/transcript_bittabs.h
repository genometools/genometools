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

#include "core/bittab.h"

/* a container class for transcript bittabs */
typedef struct TranscriptGT_Bittabs TranscriptGT_Bittabs;

/* create an empy container */
TranscriptGT_Bittabs* transcript_bittabs_new(unsigned long size_all,
                                          unsigned long size_single,
                                          unsigned long size_initial,
                                          unsigned long size_internal,
                                          unsigned long size_terminal);

/* return the bittab for all exons */
GT_Bittab*            transcript_bittabs_get_all(const TranscriptGT_Bittabs*);

/* return the bittab for single exons */
GT_Bittab*            transcript_bittabs_get_single(const TranscriptGT_Bittabs*);

/* return the bittab for initial exons */
GT_Bittab*            transcript_bittabs_get_initial(const TranscriptGT_Bittabs*);

/* return the bittab for internal exons */
GT_Bittab*            transcript_bittabs_get_internal(const TranscriptGT_Bittabs*);

/* return the bittab for terminal exons */
GT_Bittab*            transcript_bittabs_get_terminal(const TranscriptGT_Bittabs*);

void               transcript_bittabs_delete(TranscriptGT_Bittabs*);

#endif
