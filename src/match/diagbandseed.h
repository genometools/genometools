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

#ifndef SEED_EXTEND_H
#define SEED_EXTEND_H
#include <stdbool.h>
#include "core/arraydef.h"
#include "core/encseq_api.h"
#include "core/types_api.h"

struct GtSeedExtend {
  unsigned int kmerlen;
  unsigned int diagbandw;
  unsigned int mincoverage;
  unsigned int maxfreq;
  unsigned int correlation;
  unsigned int minalilen;
  unsigned int maxfrontdist;
  unsigned int minquality;
  bool mirror;
  bool verify;
  bool benchmark;
};

typedef struct GtSeedExtend GtSeedExtend;
typedef struct GtSeedExtendKmerPos GtSeedExtendKmerPos;
typedef struct GtSeedExtendSeedPair GtSeedExtendSeedPair;
GT_DECLAREARRAYSTRUCT(GtSeedExtendSeedPair);

/* Returns a GtSeedExtendKmerPos list of k-mers from a given encseq. */
GtUword gt_seed_extend_get_kmers(GtSeedExtendKmerPos *list,
                                 const GtEncseq *encseq, unsigned int k);

/* Returns a GtSeedExtendSeedPair list of equal k-mers from lists a and b. */
void gt_seed_extend_merge(GtArrayGtSeedExtendSeedPair *mlist,
                          const GtSeedExtendKmerPos *alist, GtUword alen,
                          const GtSeedExtendKmerPos *blist, GtUword blen,
                          unsigned int kmerlen, unsigned int maxfreq);

/* reports seeds from mlist that satisfy the filter criteria */
void gt_seed_extend_find_seeds(const GtArrayGtSeedExtendSeedPair *mlist,
                               unsigned int kmerlen,
                               unsigned int mincoverage,
                               unsigned int diagbandw,
                               GtUword amaxlen,
                               GtUword bmaxlen);

/* Run the whole algorithm. */
void gt_seed_extend_run(const GtEncseq *aencseq, const GtEncseq *bencseq,
                        const GtSeedExtend *arg);
#endif
