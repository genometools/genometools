/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef DP_OPTIONS_POSTPRO_H
#define DP_OPTIONS_POSTPRO_H

#include "gth/cutoffmode.h"

/* These paramters are used for post processing the ``raw'' spliced alignments
   after the DP. */
typedef struct {
  GthCutoffmode leadcutoffsmode,         /* leading cutoffs mode */
                termcutoffsmode;         /* terminal cutoffs mode */
  unsigned int cutoffsminexonlen,     /* minimum exon length for cutoffs
                                         determination */
               scoreminexonlen;       /* minimum exon length for score
                                         determination */
} GthDPOptionsPostpro;

GthDPOptionsPostpro* gth_dp_options_postpro_new(void);
void                 gth_dp_options_postpro_delete(GthDPOptionsPostpro*);

#endif
