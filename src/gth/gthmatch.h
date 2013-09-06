/*
  Copyright (c) 2004-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef GTHMATCH_H
#define GTHMATCH_H

#include <stdbool.h>
#include "core/types_api.h"

/*
  The following structure is used to store matches in gth which have been
  computed by vmatch. It is a subset of the StoreMatch structure.
  Some parts have been deleted to reduce the space consumption of gth.

  Keep in sync with function gth_matches_are_equal() in file gthmatch.c!
*/

typedef struct {
  GtWord        Storescore;             /* the score */
  GtUword       Storepositionreference, /* position of reference instance
                                           w.r.t. sequence start */
                Storelengthreference,   /* length of reference instance */
                Storepositiongenomic,   /* position of genomic instance
                                           w.r.t. sequence start */
                Storelengthgenomic,     /* length of genomic match instance */
                Storeseqnumreference,   /* sequence number of reference match
                                           instance */
                Storeseqnumgenomic;     /* sequence number of genomic match
                                           instance */
} GthMatch;

bool gth_matches_are_equal(GthMatch *match1, GthMatch *match2);

#endif
