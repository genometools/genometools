/*
  Copyright (c) 2015 Joerg Winkler <joerg.winkler@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#ifndef DIAGBANDSEED_H
#define DIAGBANDSEED_H
#include <stdbool.h>
#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/types_api.h"
#include "match/seed-extend.h"

typedef struct {
  GtUword errorpercentage,
          userdefinedleastlength;
  unsigned int seedlength;
  GtUword logdiagbandwidth;
  GtUword mincoverage;
  GtUword maxfreq;
  GtUword memlimit;
  bool norev;
  bool nofwd;
  bool overlappingseeds;
  bool verify;
  bool verbose;
  bool debug_kmer;
  bool debug_seedpair;
  bool seed_display;
  bool extend_last;
  GtGreedyextendmatchinfo *extendgreedyinfo;
  GtXdropmatchinfo *extendxdropinfo;
  GtQuerymatchoutoptions *querymatchoutopt;
} GtDiagbandseed;

/* Run the whole algorithm. */
int gt_diagbandseed_run(const GtEncseq *aencseq,
                        const GtEncseq *bencseq,
                        const GtDiagbandseed *arg,
                        GtError *err);
#endif
