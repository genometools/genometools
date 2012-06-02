/*
  Copyright (c) 2011-2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_CONTFINDER_H
#define RDJ_CONTFINDER_H

#include <stdint.h>
#include "core/file.h"
#include "match/reads2twobit.h"

typedef struct GtContfinder GtContfinder;

GtContfinder* gt_contfinder_new(GtReads2Twobit *r2t);

void gt_contfinder_run(GtContfinder *contfinder, bool mirrored,
    bool calculate_copynum);

GtBitsequence *gt_contfinder_contained(GtContfinder *contfinder);

unsigned long gt_contfinder_nofcontained(GtContfinder *contfinder);

int gt_contfinder_write_seqnums(GtContfinder *contfinder, bool sorted,
    GtFile *outfp, GtError *err);

int gt_contfinder_write_sorted_seqnums(GtContfinder *contfinder, char* path,
    GtError *err);

int gt_contfinder_write_cntlist(GtContfinder *contfinder, char* path,
    GtError *err);

int gt_contfinder_write_copynum(GtContfinder *contfinder, char* path,
    GtError *err);

void gt_contfinder_delete(GtContfinder *contfinder);

void gt_contfinder_radixsort_str_eqlen_tester(GtContfinder *contfinder,
    bool mirrored, unsigned long depth,
    unsigned long maxdepth, bool print);

#endif
