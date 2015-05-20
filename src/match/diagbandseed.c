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
#include "core/unused_api.h"
#include "core/alphabet_api.h"

#define GT_SEED_EXTEND_ARRAY_INCR 256

typedef uint32_t GtSeedExtendPosition;
typedef uint32_t GtSeedExtendSeqnum;

struct GtSeedExtendKmerPos {
  GtCodetype code;            /* only sort criterion */
  GtSeedExtendSeqnum seqnum;
  GtSeedExtendPosition endpos;
};

struct GtSeedExtendSeedPair {
  GtSeedExtendSeqnum bseqnum; /*  2nd important sort criterion */
  GtSeedExtendSeqnum aseqnum; /* most important sort criterion */
  GtSeedExtendPosition apos;
  GtSeedExtendPosition bpos;  /*  3rd important sort criterion */
};

/* Returns a GtSeedExtendKmerPos list of k-mers from a given encseq. */
void gt_seed_extend_get_kmers(GtArrayGtSeedExtendKmerPos *list,
                              const GtEncseq *encseq, unsigned int kmerlen)
{
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;

  gt_assert(list != NULL && encseq != NULL);
  kc_iter = gt_kmercodeiterator_encseq_new(encseq, GT_READMODE_FORWARD, kmerlen,
                                           0);
  while ((kmercode = gt_kmercodeiterator_encseq_next(kc_iter)) != NULL) {
    if (!kmercode->definedspecialposition) {
      /* store (kmercode, seqnum, endpos) in array */
      GtSeedExtendKmerPos *kmerposptr;
      GT_GETNEXTFREEINARRAY(kmerposptr, list, GtSeedExtendKmerPos,
                            GT_SEED_EXTEND_ARRAY_INCR);
      kmerposptr->code = kmercode->code;
      kmerposptr->endpos = gt_kmercodeiterator_encseq_get_currentpos(kc_iter)-1;
      kmerposptr->seqnum = gt_encseq_seqnum(encseq, kmerposptr->endpos);
      gt_assert(kmerposptr->endpos >= gt_encseq_seqstartpos(encseq,
                kmerposptr->seqnum));
      kmerposptr->endpos -= gt_encseq_seqstartpos(encseq, kmerposptr->seqnum);
#ifdef MYTEST
      if (true) {
        char *buf = gt_malloc(kmerlen*sizeof(char));
        GtSeedExtendPosition tmp
          = gt_kmercodeiterator_encseq_get_currentpos(kc_iter) - 1;
        gt_encseq_extract_decoded(encseq, buf, tmp+1-kmerlen, tmp);
        printf("Kmer (%7lu, %s, %d, %d)\n", kmercode->code, buf,
               kmerposptr->seqnum, kmerposptr->endpos);
        gt_free(buf);
      }
#endif
    } else {
      /* if specialposition is N: store kmer in array */
      /* allow comparison of N-containing kmers */
    }
  }
  gt_kmercodeiterator_delete(kc_iter);
}

/* Returns a GtSeedExtendSeedPair list of equal kmers from lists a and b. */
void gt_seed_extend_merge(GtArrayGtSeedExtendSeedPair *mlist,
                          const GtArrayGtSeedExtendKmerPos *alist,
                          const GtArrayGtSeedExtendKmerPos *blist)
{
  const GtSeedExtendKmerPos *aptr, *bptr, *tptr, *aend, *bend;
  gt_assert(alist != NULL && blist != NULL && mlist != NULL);
  aptr = alist->spaceGtSeedExtendKmerPos;
  bptr = blist->spaceGtSeedExtendKmerPos;
  aend = aptr + alist->nextfreeGtSeedExtendKmerPos;
  bend = bptr + blist->nextfreeGtSeedExtendKmerPos;
  while (aptr < aend && bptr < bend) {
    if (aptr->code < bptr->code) {
      aptr ++;
    } else if (aptr->code > bptr->code) {
      bptr ++;
    } else { /* k-mer codes are equal: process all equal elements from blist */
      for (tptr = bptr; tptr < bend && aptr->code == tptr->code; tptr ++) {
        if (alist != blist || aptr->seqnum < tptr->seqnum) { /* no duplicates */
          GtSeedExtendSeedPair *seedptr;
          GT_GETNEXTFREEINARRAY(seedptr, mlist, GtSeedExtendSeedPair,
                                GT_SEED_EXTEND_ARRAY_INCR);
          seedptr->bseqnum = tptr->seqnum;
          seedptr->aseqnum = aptr->seqnum;
          seedptr->bpos = tptr->endpos;
          seedptr->apos = aptr->endpos;
        }
      }
      aptr ++;
    }
  }
}

bool gt_seed_extend_is_seed(GT_UNUSED const GtSeedExtendSeedPair *entry,
                            const unsigned int *score, unsigned int mincoverage,
                            unsigned int diag)
{
  /* TODO: add path criteron */
  gt_assert(score != NULL);
  if (MAX(score[diag+1], score[diag-1]) + score[diag] >= mincoverage) {
    return true;
  } else {
    return false;
  }
}

void gt_seed_extend_find_seeds(const GtArrayGtSeedExtendSeedPair *mlist,
                               unsigned int kmerlen, unsigned int mincoverage,
                               unsigned int diagbandw, GtUword amaxlen,
                               GtUword bmaxlen)
{
  const GtUword mlen = mlist->nextfreeGtSeedExtendSeedPair;
  const GtSeedExtendSeedPair *lm = mlist->spaceGtSeedExtendSeedPair;
  const GtUword ndiags = (amaxlen >> diagbandw) + (bmaxlen >> diagbandw) + 2;
  const GtUword minhit = (mincoverage-1) / kmerlen + 1;
  unsigned int *score = gt_calloc(ndiags, sizeof *score);
  GtSeedExtendPosition *lastp = gt_calloc(ndiags, sizeof *lastp);
  GtUword diag, nextsegm, i;

#ifdef MYTEST
  printf("ndiags="GT_WU", s=%d\n", ndiags, diagbandw);
  printf("amaxlen="GT_WU", bmaxlen="GT_WU"\n", amaxlen, bmaxlen);
#endif

  nextsegm = 0;
  while (nextsegm + minhit <= mlen) {
    const GtUword currsegm = nextsegm;
    /* if insuffienct number of kmers in segment: skip whole segment */
    if (lm[nextsegm].aseqnum != lm[nextsegm+minhit-1].aseqnum ||
        lm[nextsegm].bseqnum != lm[nextsegm+minhit-1].bseqnum) {
      do {
        nextsegm ++;
      } while (nextsegm < mlen &&
               lm[nextsegm].aseqnum == lm[currsegm].aseqnum &&
               lm[nextsegm].bseqnum == lm[currsegm].bseqnum);
      continue;
    }

    /* calculate diagonal band scores */
    do {
      gt_assert(lm[nextsegm].bpos <= bmaxlen && lm[nextsegm].apos <= amaxlen);
      diag = (amaxlen + lm[nextsegm].bpos - lm[nextsegm].apos) >> diagbandw;
      if (lm[nextsegm].bpos >= kmerlen + lastp[diag]) {
        score[diag] += kmerlen;
      } else {
        gt_assert(lastp[diag] <= lm[nextsegm].bpos);/*if fail: sorted by bpos?*/
        score[diag] = score[diag] + lm[nextsegm].bpos - lastp[diag];
      }
      lastp[diag] = lm[nextsegm].bpos;
      nextsegm ++;
    } while (nextsegm < mlen &&
             lm[nextsegm].aseqnum == lm[currsegm].aseqnum &&
             lm[nextsegm].bseqnum == lm[currsegm].bseqnum);

    /* report seeds */
    for (i = currsegm; i < nextsegm; i++) {
      gt_assert(lm[i].apos <= amaxlen);
      diag = (amaxlen + lm[i].bpos - lm[i].apos) >> diagbandw;
      if (gt_seed_extend_is_seed(&lm[i], score, mincoverage, diag)) {
#ifdef MYTEST
        printf("report SeedPair (%d,%d,%d,%d), score["GT_WU"]=%d\n",
               lm[i].aseqnum, lm[i].bseqnum, lm[i].apos, lm[i].bpos, diag,
               MAX(score[diag+1], score[diag-1]) + score[diag]);
#endif
      }
    }

    /* reset diagonal band scores */
    for (i = currsegm; i < nextsegm; i++) {
      diag = (amaxlen + lm[i].bpos - lm[i].apos) >> diagbandw;
      score[diag] = 0;
      lastp[diag] = 0;
    }
  }
  gt_free(score);
  gt_free(lastp);
}

void gt_seed_extend_run(const GtEncseq *aencseq, const GtEncseq *bencseq,
                        unsigned int kmerlen,
                        unsigned int mincoverage,
                        unsigned int diagbandw)
{
  GtArrayGtSeedExtendKmerPos alist, blist;
  GtArrayGtSeedExtendSeedPair mlist;
  GtRadixsortinfo* rdxinfo;
  const bool two_files = (bencseq != aencseq) ? true : false;
  const GtUword amaxlen = gt_encseq_max_seq_length(aencseq);
  const GtUword bmaxlen = gt_encseq_max_seq_length(bencseq);

  if (amaxlen < kmerlen || bmaxlen < kmerlen) {
    /* printf("maximum sequence length too short to find any kmers\n"); */
    return;
  }

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
    gt_radixsort_inplace_GtUwordPair((GtUwordPair*)blist.
                                     spaceGtSeedExtendKmerPos,
                                     blist.nextfreeGtSeedExtendKmerPos);
    gt_radixsort_delete(rdxinfo);
  } else { /* compare all reads of encseq A with themselves */
    /* bencseq = aencseq; */
  }

#ifdef MYTEST
  if (alist.nextfreeGtSeedExtendKmerPos != 0) {
    char *buf = gt_malloc(kmerlen*sizeof(char));
    GtSeedExtendKmerPos *j = alist.spaceGtSeedExtendKmerPos;
    GtSeedExtendKmerPos *last = j + alist.nextfreeGtSeedExtendKmerPos;
    while (j < last) {
      gt_encseq_extract_decoded(aencseq, buf, j->endpos+1-kmerlen, j->endpos);
      printf("Kmer (%lu, %d, %d)\n", j->code, j->seqnum, j->endpos);
      j ++;
    }
    gt_free(buf);
  }
#endif

  GT_INITARRAY(&mlist,GtSeedExtendSeedPair);
  if (two_files)
    gt_seed_extend_merge(&mlist, &alist, &blist);
  else
    gt_seed_extend_merge(&mlist, &alist, &alist);
  rdxinfo = gt_radixsort_new_uint64keypair(mlist.nextfreeGtSeedExtendSeedPair);
  gt_radixsort_inplace_Gtuint64keyPair((Gtuint64keyPair*)mlist.
                                       spaceGtSeedExtendSeedPair,
                                       mlist.nextfreeGtSeedExtendSeedPair);
  gt_radixsort_delete(rdxinfo);

#ifndef MYTEST
  if (mlist.nextfreeGtSeedExtendSeedPair != 0) {
    GtSeedExtendSeedPair *j = mlist.spaceGtSeedExtendSeedPair;
    GtSeedExtendSeedPair *last = j + mlist.nextfreeGtSeedExtendSeedPair;
    char *buf = gt_malloc(kmerlen*sizeof(char));
    char *buf2 = gt_malloc(kmerlen*sizeof(char));
    while (j < last) {
      GtSeedExtendPosition a=j->apos+gt_encseq_seqstartpos(aencseq, j->aseqnum);
      GtSeedExtendPosition b=j->bpos+gt_encseq_seqstartpos(bencseq, j->bseqnum);
      gt_encseq_extract_decoded(aencseq, buf, a+1-kmerlen, a);
      gt_encseq_extract_decoded(bencseq, buf2, b+1-kmerlen, b);
      printf("SeedPair (%d,%d,%d,%d)\n", j->aseqnum, j->bseqnum, j->apos,
             j->bpos);
      j ++;
    }
    gt_free(buf);
    gt_free(buf2);
  }
#endif

  if (mlist.nextfreeGtSeedExtendSeedPair != 0)
    gt_seed_extend_find_seeds(&mlist, kmerlen, mincoverage, diagbandw, amaxlen,
                              bmaxlen);

  GT_FREEARRAY(&mlist, GtSeedExtendSeedPair);
  GT_FREEARRAY(&alist, GtSeedExtendKmerPos);
  if (two_files)
    GT_FREEARRAY(&blist, GtSeedExtendKmerPos);
}
