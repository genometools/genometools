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

#include "core/codetype.h"
#include "core/encseq_api.h"
#include "core/radix_sort.h"
#include "match/seed-extend.h"
#include "match/sfx-mappedstr.h"

#define GT_SEED_EXTEND_ARRAY_INCR 256

typedef uint32_t GtSeedExtendPosition;
typedef uint32_t GtSeedExtendRead;

struct GtSeedExtendKmerPos {
  GtCodetype code;            /* only sort criterion */
  GtSeedExtendRead read;
  GtSeedExtendPosition endpos;
};

struct GtSeedExtendSeedPair {
  GtSeedExtendRead bread;     /*  2nd important sort criterion */
  GtSeedExtendRead aread;     /* most important sort criterion */
  GtSeedExtendPosition bpos;
  GtSeedExtendPosition apos;  /*  3rd important sort criterion */
};

/* Returns a GtSeedExtendKmerPos list of k-mers from a given encseq. */
void gt_seed_extend_get_kmers(GtArrayGtSeedExtendKmerPos *list,
                              const GtEncseq *encseq, const unsigned int k)
{
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;
  
  kc_iter = gt_kmercodeiterator_encseq_new(encseq, GT_READMODE_FORWARD, k, 0);
  while (kmercode = gt_kmercodeiterator_encseq_nonspecial_next(kc_iter)) {
    GtSeedExtendKmerPos kp;
    kp.code = kmercode->code;
    kp.endpos = gt_kmercodeiterator_encseq_get_currentpos(kc_iter) - 1;
    kp.read = gt_encseq_seqnum(encseq, kp.endpos);
    GT_STOREINARRAY(list, GtSeedExtendKmerPos, GT_SEED_EXTEND_ARRAY_INCR, kp);
  }
  gt_kmercodeiterator_delete(kc_iter);
}

/* Returns a GtSeedExtendSeedPair list of equal kmers from lists a and b. */
void gt_seed_extend_merge(GtArrayGtSeedExtendSeedPair *mlist,
                          const GtArrayGtSeedExtendKmerPos *alist,
                          const GtArrayGtSeedExtendKmerPos *blist)
{
  const GtSeedExtendKmerPos *aptr, *bptr, *tptr, *aend, *bend;
  aptr = alist->spaceGtSeedExtendKmerPos;
  bptr = blist->spaceGtSeedExtendKmerPos;
  aend = aptr + alist->nextfreeGtSeedExtendKmerPos;
  bend = bptr + blist->nextfreeGtSeedExtendKmerPos;
  while (aptr < aend && bptr < bend) {
    if (aptr->code < bptr->code) {
      aptr ++;
    } else if (aptr->code > bptr->code) {
      bptr ++;
    } else { /* k-mer codes are equal */
      for (tptr = bptr; tptr < bend && aptr->code == tptr->code; tptr ++) {
        if (alist != blist || aptr->read != tptr->read) {
          GtSeedExtendSeedPair seed;
          seed.bread = tptr->read;
          seed.aread = aptr->read;
          seed.bpos = tptr->endpos;
          seed.apos = aptr->endpos;
          GT_STOREINARRAY(mlist, GtSeedExtendSeedPair, 
                          GT_SEED_EXTEND_ARRAY_INCR, seed);
        }
      }
      aptr ++;
    }
  }
}

void gt_seed_extend_run(GtEncseq *aencseq, GtEncseq *bencseq,
                        const unsigned int kmerlen,
                        GT_UNUSED const unsigned int mincoverage,
                        GT_UNUSED const unsigned int diagbandw)
{
  GtArrayGtSeedExtendKmerPos alist, blist;
  GtArrayGtSeedExtendSeedPair mlist;
  GtRadixsortinfo* rdxinfo;
  const bool two_files = (bencseq != aencseq);

  GT_INITARRAY(&alist,GtSeedExtendKmerPos);
  gt_seed_extend_get_kmers(&alist, aencseq, kmerlen);
  rdxinfo = gt_radixsort_new_ulongpair(alist.nextfreeGtSeedExtendKmerPos);
  gt_radixsort_inplace_GtUwordPair((GtUwordPair*)alist.spaceGtSeedExtendKmerPos,
                                   alist.nextfreeGtSeedExtendKmerPos);
  gt_radixsort_delete(rdxinfo);
  
  if (two_files) {
    GT_INITARRAY(&blist,GtSeedExtendKmerPos);
    gt_seed_extend_get_kmers(&blist, bencseq, kmerlen);
    rdxinfo = gt_radixsort_new_ulongpair(blist.nextfreeGtSeedExtendKmerPos);
    gt_radixsort_inplace_GtUwordPair((GtUwordPair*)blist.spaceGtSeedExtendKmerPos,
                                     blist.nextfreeGtSeedExtendKmerPos);
    gt_radixsort_delete(rdxinfo);
  } else { /* compare all reads of encseq A with themselves */
    bencseq = aencseq;
  }
  
  GT_INITARRAY(&mlist,GtSeedExtendSeedPair);
  if (two_files)
    gt_seed_extend_merge(&mlist, &alist, &blist);
  else
    gt_seed_extend_merge(&mlist, &alist, &alist);
  rdxinfo = gt_radixsort_new_ulongpair(mlist.nextfreeGtSeedExtendSeedPair);
  gt_radixsort_inplace_GtUwordPair((GtUwordPair*)mlist.spaceGtSeedExtendSeedPair,
                                   mlist.nextfreeGtSeedExtendSeedPair);
  gt_radixsort_delete(rdxinfo);
  
#ifndef NOPRINT
  if (mlist.nextfreeGtSeedExtendSeedPair != 0) {
    GtSeedExtendSeedPair *j = mlist.spaceGtSeedExtendSeedPair;
    GtSeedExtendSeedPair *last = j + mlist.nextfreeGtSeedExtendSeedPair;
    while (j < last) {
      printf("SeedPair (%d,%d,%d,%d)\n", j->aread, j->bread, j->apos, j->bpos);
      j ++;
    }
  }
#endif

  GT_FREEARRAY(&mlist, GtSeedExtendSeedPair);
  GT_FREEARRAY(&alist, GtSeedExtendKmerPos);
  if (two_files)
    GT_FREEARRAY(&blist, GtSeedExtendKmerPos);
}

