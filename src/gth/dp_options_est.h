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

#ifndef DP_OPTIONS_EST_H
#define DP_OPTIONS_EST_H

#include <stdbool.h>

typedef struct {
  /* basic DP options */
  double probies;           /* P that initial state is exon state */
  double probdelgen,        /* P of deletion in genomic sequence */
         identityweight,    /* identity (sigma) */
         mismatchweight,    /* mismatch (mu) */
         undetcharweight,   /* alignment positions involving undetermined
                               characters (nu) */
         deletionweight;    /* cDNA deletions (delta) */

  /* special options */
  unsigned int wzerotransition,  /* window size for zero transition weights */
               wdecreasedoutput; /* window size for decreased output weights */

  /* development options */
  bool detectsmallexons;         /* try to detect small exons */
} GthDPOptionsEST;

GthDPOptionsEST* gth_dp_options_est_new(void);
GthDPOptionsEST* gth_dp_options_est_clone(const GthDPOptionsEST*);
void             gth_dp_options_est_delete(GthDPOptionsEST*);

#endif
