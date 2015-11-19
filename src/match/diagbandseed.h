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
#include "core/arraydef.h"
#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/readmode_api.h"
#include "core/types_api.h"
#include "match/seed-extend.h"
#include "match/xdrop.h"

struct GtDiagbandseed {
  GtUword errorpercentage,
          userdefinedleastlength;
  unsigned int seedlength;
  GtUword logdiagbandwidth;
  GtUword mincoverage;
  GtUword maxfreq;
  GtUword memlimit;
  bool mirror;
  bool overlappingseeds;
  bool verify;
  bool verbose;
  bool debug_kmer;
  bool debug_seedpair;
  bool seed_display;
  GtGreedyextendmatchinfo *extendgreedyinfo;
  GtXdropmatchinfo *extendxdropinfo;
  GtQuerymatchoutoptions *querymatchoutopt;
};

typedef struct GtDiagbandseed GtDiagbandseed;
typedef struct GtDiagbandseedKmerPos GtDiagbandseedKmerPos;
typedef struct GtDiagbandseedSeedPair GtDiagbandseedSeedPair;
GT_DECLAREARRAYSTRUCT(GtDiagbandseedSeedPair);

/* Returns a GtDiagbandseedKmerPos list of k-mers from a given encseq. */
int gt_diagbandseed_get_kmers(GtDiagbandseedKmerPos *list,
                              GtUword *listlength,
                              const GtEncseq *encseq,
                              unsigned int seedlength,
                              GtReadmode readmode,
                              GtError *err);

/* Returns a GtDiagbandseedSeedPair list of equal k-mers from lists a and b. */
void gt_diagbandseed_merge(GtArrayGtDiagbandseedSeedPair *mlist,
                           const GtDiagbandseedKmerPos *alist, GtUword alen,
                           const GtDiagbandseedKmerPos *blist, GtUword blen,
                           GtUword *maxfreq,
                           GtUword maxgram,
                           GtUword memlimit,
                           GtUword *histogram,
                           unsigned int endposdiff,
                           bool selfcomp);

/* Run the whole algorithm. */
int gt_diagbandseed_run(const GtEncseq *aencseq,
                        const GtEncseq *bencseq,
                        const GtDiagbandseed *arg,
                        GtError *err);
#endif
