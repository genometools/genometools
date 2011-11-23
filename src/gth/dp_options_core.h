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

#ifndef DP_OPTIONS_CORE_H
#define DP_OPTIONS_CORE_H

#include <stdbool.h>
#include "core/range_api.h"

typedef struct {
  bool noicinintroncheck,         /* perform no check if intron coutout is in
                                     intron */
       freeintrontrans;           /* free state transitions between intron
                                     states */
  unsigned int dpminexonlength,   /* minimum exon length for the DP */
               dpminintronlength; /* minimum intron length */
  double shortexonpenalty,        /* penalty for short exons */
         shortintronpenalty;      /* penalty for short introns */
  GtRange btmatrixgenrange,
          btmatrixrefrange;
  unsigned long jtoverlap;
  bool jtdebug;
} GthDPOptionsCore;

GthDPOptionsCore* gth_dp_options_core_new(void);
GthDPOptionsCore* gth_dp_options_core_clone(const GthDPOptionsCore*);
void              gth_dp_options_core_delete(GthDPOptionsCore *);

#endif
