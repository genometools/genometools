/*
 Copyright (c) 2015 Jörg Winkler <joerg.winkler@studium.uni-hamburg.de>
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
#include "core/types_api.h"

struct GtDiagbandseed {
  unsigned int dbs_seedlength;
  GtUword dbs_logdiagbandwidth;
  GtUword dbs_mincoverage;
  GtUword dbs_maxfreq;
  GtUword se_alignlength;
  GtUword se_maxalilendiff;
  GtUword se_errorpercentage;
  GtUword se_historysize;
  bool mirror;
  bool overlappingseeds;
  bool verify;
  bool benchmark;
};

typedef struct GtDiagbandseed GtDiagbandseed;
typedef struct GtDiagbandseedKmerPos GtDiagbandseedKmerPos;
typedef struct GtDiagbandseedSeedPair GtDiagbandseedSeedPair;
GT_DECLAREARRAYSTRUCT(GtDiagbandseedSeedPair);

/* Returns a GtDiagbandseedKmerPos list of k-mers from a given encseq. */
GtUword gt_diagbandseed_get_kmers(GtDiagbandseedKmerPos *list,
                                  const GtEncseq *encseq, unsigned int kmerlen,
                                  bool revcompl);

/* Returns a GtDiagbandseedSeedPair list of equal k-mers from lists a and b. */
void gt_diagbandseed_merge(GtArrayGtDiagbandseedSeedPair *mlist,
                           const GtDiagbandseedKmerPos *alist, GtUword alen,
                           const GtDiagbandseedKmerPos *blist, GtUword blen,
                           unsigned int kmerlen, GtUword maxfreq,
                           bool two_files);

/* reports seeds from mlist that satisfy the filter criteria */
void gt_diagbandseed_find_seeds(const GtArrayGtDiagbandseedSeedPair *mlist,
                                unsigned int kmerlen, GtUword mincoverage,
                                GtUword log_diagbandwidth, GtUword amaxlen,
                                GtUword bmaxlen);

/* Run the whole algorithm. */
void gt_diagbandseed_run(const GtEncseq *aencseq, const GtEncseq *bencseq,
                         const GtDiagbandseed *arg);
#endif
