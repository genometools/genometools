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
#include "core/minmax.h"
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
    kp.endpos -= gt_encseq_seqstartpos(encseq, gt_encseq_seqnum(encseq, 
                                                                kp.endpos));
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

bool gt_seed_extend_is_seed(GT_UNUSED const GtSeedExtendSeedPair *entry, 
                            const unsigned int *score,
                            const unsigned int mincoverage,
                            const unsigned int diag)
{
  // TODO: add path criteron
  if (MAX(score[diag+1], score[diag-1]) + score[diag] >= mincoverage) {
    return true;
  } else {
    return false;
  }
}

void gt_seed_extend_find_seeds(const GtArrayGtSeedExtendSeedPair *mlist,
                               const unsigned int kmerlen,
                               const unsigned int mincoverage,
                               const unsigned int diagbandw,
                               const GtUword amaxlen, const GtUword bmaxlen)
{
  const GtUword mlen = mlist->nextfreeGtSeedExtendSeedPair;
  const GtSeedExtendSeedPair *lm = mlist->spaceGtSeedExtendSeedPair;
  const GtUword offset = (amaxlen >> diagbandw) + 1;
  const GtUword ndiags = (bmaxlen >> diagbandw) + 1 + offset;
  const GtUword minhit = (mincoverage-1) / kmerlen + 1;
  unsigned int *score = gt_calloc(ndiags, sizeof(unsigned int));
  GtSeedExtendPosition *lastp = gt_calloc(ndiags, sizeof(GtSeedExtendPosition));
  GtUword diag;
  unsigned int i, f;
  
#ifdef MYTEST
  printf("ndiags="GT_WU", offset="GT_WU", s=%d\n", ndiags, offset, diagbandw);
  printf("amaxlen="GT_WU", bmaxlen="GT_WU"\n", amaxlen, bmaxlen);
#endif
  
  i = 0;
  while (i+minhit <= mlen) {
    const unsigned int p = i;
    /* if insuffienct number of kmers in segment: skip whole segment */
    if (lm[i].aread != lm[i+minhit-1].aread || 
        lm[i].bread != lm[i+minhit-1].bread) {
      do {
        i ++;
      } while (i < mlen && lm[i].aread == lm[p].aread && lm[i].bread == lm[p].bread);
      continue;
    }
    
    /* calculate diagonal band scores */
    do {
      diag = offset + (((long)lm[i].apos - (long)lm[i].bpos) >> diagbandw);
      if (lm[i].apos >= kmerlen + lastp[diag]) {
        score[diag] += kmerlen;
      } else {
        score[diag] += lm[i].apos - lastp[diag];
      }
      lastp[diag] = lm[i].apos;
      i++;
    } while (i < mlen && lm[i].aread == lm[p].aread && lm[i].bread == lm[p].bread);
    
    /* report seeds */
    for (f = p; f < i; f++) {
      diag = offset + (((long)lm[f].apos - (long)lm[f].bpos) >> diagbandw);
      if (gt_seed_extend_is_seed(&lm[f], score, mincoverage, diag)) {
#ifdef MYTEST
        printf("report SeedPair (%d,%d,%d,%d), score["GT_WU"]=%d\n", 
               lm[f].aread, lm[f].bread, lm[f].apos, lm[f].bpos, diag,
               MAX(score[diag+1], score[diag-1]) + score[diag]);
#endif
      }
    }
    
    /* reset diagonal band scores */
    for (f = p; f < i; f++) {
      diag = offset + (((long)lm[f].apos - (long)lm[f].bpos) >> diagbandw);
      score[diag] = 0;
      lastp[diag] = 0;
    }
  }
  gt_free(score);
  gt_free(lastp);
}

void gt_seed_extend_run(GtEncseq *aencseq, GtEncseq *bencseq,
                        const unsigned int kmerlen,
                        const unsigned int mincoverage,
                        const unsigned int diagbandw)
{
  GtArrayGtSeedExtendKmerPos alist, blist;
  GtArrayGtSeedExtendSeedPair mlist;
  GtRadixsortinfo* rdxinfo;
  const bool two_files = (bencseq != aencseq);
  const GtUword amaxlen = gt_encseq_max_seq_length(aencseq);
  const GtUword bmaxlen = gt_encseq_max_seq_length(bencseq);

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
  
  gt_seed_extend_find_seeds(&mlist, kmerlen, mincoverage, diagbandw, amaxlen, 
                            bmaxlen);

  GT_FREEARRAY(&mlist, GtSeedExtendSeedPair);
  GT_FREEARRAY(&alist, GtSeedExtendKmerPos);
  if (two_files)
    GT_FREEARRAY(&blist, GtSeedExtendKmerPos);
}

