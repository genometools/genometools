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

#include <string.h>
#include <stdio.h>
#include "core/alphabet_api.h"
#include "core/codetype.h"
#include "core/minmax.h"
#include "core/radix_sort.h"
#include "core/timer_api.h"
#include "core/unused_api.h"
#include "match/seed-extend.h"
#include "match/sfx-mappedstr.h"

#define DEBUG_SEEDPAIR
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
GtUword gt_seed_extend_get_kmers(GtSeedExtendKmerPos *list,
                              const GtEncseq *encseq, unsigned int kmerlen)
{
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;
  GtUword length = 0;

  gt_assert(list != NULL && encseq != NULL);
  kc_iter = gt_kmercodeiterator_encseq_new(encseq, GT_READMODE_FORWARD, kmerlen,
                                           0);
  while ((kmercode = gt_kmercodeiterator_encseq_next(kc_iter)) != NULL) {
    if (!kmercode->definedspecialposition) {
      /* store (kmercode, seqnum, endpos) in array */
      GtSeedExtendKmerPos *kmerposptr = list + length;
      kmerposptr->code = kmercode->code;
      kmerposptr->endpos = gt_kmercodeiterator_encseq_get_currentpos(kc_iter)-1;
      kmerposptr->seqnum = gt_encseq_seqnum(encseq, kmerposptr->endpos);
      gt_assert(kmerposptr->endpos >= gt_encseq_seqstartpos(encseq,
                kmerposptr->seqnum));
      kmerposptr->endpos -= gt_encseq_seqstartpos(encseq, kmerposptr->seqnum);
      length++;
    } else {
      /* if specialposition is N: store kmer in array */
      /* allow comparison of N-containing kmers */
    }
  }
  gt_kmercodeiterator_delete(kc_iter);
  return length;
}

/* Returns a GtSeedExtendSeedPair list of equal kmers from lists a and b. */
void gt_seed_extend_merge(GtArrayGtSeedExtendSeedPair *mlist,
                          const GtSeedExtendKmerPos *alist, GtUword alen,
                          const GtSeedExtendKmerPos *blist, GtUword blen,
                          unsigned int maxfreq)
{
  const GtSeedExtendKmerPos *aptr, *bptr, *aend, *bend;
  gt_assert(alist != NULL && blist != NULL && mlist != NULL);

  aptr = alist;
  bptr = blist;
  aend = aptr + alen;
  bend = bptr + blen;
  while (aptr < aend && bptr < bend) {
    if (aptr->code < bptr->code) {
      aptr ++;
    } else if (aptr->code > bptr->code) {
      bptr ++;
    } else {
      /* equality: count frequency of current k-mer in both lists */
      const GtSeedExtendKmerPos *aiter, *biter;
      for (aiter = aptr; aiter < aend && aiter->code == bptr->code; aiter++) {}
      for (biter = bptr; biter < bend && biter->code == aptr->code; biter++) {}
      if (aiter - aptr <= maxfreq && biter - bptr <= maxfreq) {
        /* add all equal k-mers */
        const GtSeedExtendKmerPos *asegm_end = aiter, *bsegm_end = biter;
        for (aiter = aptr; aiter < asegm_end; aiter++) {
          for (biter = bptr; biter < bsegm_end; biter++) {
            if (alist != blist || aiter->seqnum < biter->seqnum) {
              /* no duplicates from the same dataset */
              GtSeedExtendSeedPair *seedptr;
              GT_GETNEXTFREEINARRAY(seedptr, mlist, GtSeedExtendSeedPair,
                                    GT_SEED_EXTEND_ARRAY_INCR + 1.2 *
                                    mlist->allocatedGtSeedExtendSeedPair);
              seedptr->bseqnum = biter->seqnum;
              seedptr->aseqnum = aiter->seqnum;
              seedptr->bpos = biter->endpos;
              seedptr->apos = aiter->endpos;
            }
          }
        }
      } /* else: ignore all equal elements */
      aptr = aiter;
      bptr = biter;
    }
  }
}

bool gt_seed_extend_is_seed(GT_UNUSED const GtSeedExtendSeedPair *entry,
                            const unsigned int *score, unsigned int mincoverage,
                            unsigned int diag)
{
  /* TODO: add path criterion */
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
#ifdef DEBUG_SEED_REPORT
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
                        unsigned int kmerlen, unsigned int mincoverage,
                        unsigned int diagbandw, unsigned int maxfreq,
                        bool verify, bool benchmark)
{
  GtSeedExtendKmerPos *alist, *blist;
  GtArrayGtSeedExtendSeedPair mlist;
  GtRadixsortinfo* rdxinfo;
  GtUword alen, blen;
  const bool two_files = (bencseq != aencseq) ? true : false;
  const bool one_seq = (!two_files && gt_encseq_num_of_sequences(aencseq) <= 1);
  const GtUword amaxlen = gt_encseq_max_seq_length(aencseq);
  const GtUword bmaxlen = gt_encseq_max_seq_length(bencseq);
  const GtUword ankmers = gt_encseq_total_length(aencseq) -
                          MIN(kmerlen-1, gt_encseq_min_seq_length(aencseq)) *
                          gt_encseq_num_of_sequences(aencseq);
  const GtUword bnkmers = gt_encseq_total_length(bencseq) -
                          MIN(kmerlen-1, gt_encseq_min_seq_length(bencseq)) *
                          gt_encseq_num_of_sequences(bencseq);
  GtTimer *timer;

  if (one_seq || amaxlen < kmerlen || bmaxlen < kmerlen) {
    /*printf("maximum sequence length too short or only 1 sequence given\n");*/
    if (benchmark) {
      printf("0.000000,0,0\n");
    }
    return;
  }

  alist = gt_malloc(ankmers * sizeof *alist);
  alen = gt_seed_extend_get_kmers(alist, aencseq, kmerlen);
  if (benchmark) {
    timer = gt_timer_new();
    gt_timer_start(timer);
  }
  rdxinfo = gt_radixsort_new_ulongpair(alen);
  gt_radixsort_inplace_GtUwordPair((GtUwordPair*)alist, alen);
  gt_radixsort_delete(rdxinfo);

  if (two_files) {
    blist = gt_malloc(bnkmers * sizeof *blist);
    blen = gt_seed_extend_get_kmers(blist, bencseq, kmerlen);
    rdxinfo = gt_radixsort_new_ulongpair(blen);
    gt_radixsort_inplace_GtUwordPair((GtUwordPair*)blist, blen);
    gt_radixsort_delete(rdxinfo);
  } else {
    /* compare all reads of encseq A with themselves */
    blist = alist;
    blen = alen;
  }

  GT_INITARRAY(&mlist,GtSeedExtendSeedPair);
  gt_seed_extend_merge(&mlist, alist, alen, blist, blen, maxfreq);
  gt_free(alist);
  if (two_files)
    gt_free(blist);
  rdxinfo = gt_radixsort_new_uint64keypair(mlist.nextfreeGtSeedExtendSeedPair);
  gt_radixsort_inplace_Gtuint64keyPair((Gtuint64keyPair*)mlist.
                                       spaceGtSeedExtendSeedPair,
                                       mlist.nextfreeGtSeedExtendSeedPair);
  gt_radixsort_delete(rdxinfo);
  if (benchmark) {
    gt_timer_stop(timer);
    gt_timer_show_formatted(timer, GT_WD ".%06ld,"GT_WD","GT_WD"\n", stdout);
  }

  if (verify && mlist.nextfreeGtSeedExtendSeedPair != 0) {
    GtSeedExtendSeedPair *j = mlist.spaceGtSeedExtendSeedPair;
    GtSeedExtendSeedPair *last = j + mlist.nextfreeGtSeedExtendSeedPair;
    char *buf1 = gt_malloc(1 + kmerlen * sizeof *buf1);
    char *buf2 = gt_malloc(1 + kmerlen * sizeof *buf2);
    while (j < last) {
      GtSeedExtendPosition a=j->apos+gt_encseq_seqstartpos(aencseq, j->aseqnum);
      GtSeedExtendPosition b=j->bpos+gt_encseq_seqstartpos(bencseq, j->bseqnum);
      gt_encseq_extract_decoded(aencseq, buf1, a + 1 - kmerlen, a);
      gt_encseq_extract_decoded(bencseq, buf2, b + 1 - kmerlen, b);
      buf1[kmerlen] = buf2[kmerlen] = '\0';
#ifdef DEBUG_SEEDPAIR
      printf("SeedPair (%d,%d,%d,%d)\n", j->aseqnum, j->bseqnum, j->apos,
             j->bpos);
#endif
      if (strcmp(buf1, buf2) != 0) {
        fprintf(stderr, "wrong seed(%d,%d,%d,%d): %s != %s\n",
                j->aseqnum, j->bseqnum, j->apos, j->bpos, buf1, buf2);
        gt_assert(strcmp(buf1, buf2) == 0);
      }
      j++;
    }
    gt_free(buf1);
    gt_free(buf2);
  }

  if (mlist.nextfreeGtSeedExtendSeedPair != 0)
    gt_seed_extend_find_seeds(&mlist, kmerlen, mincoverage, diagbandw, amaxlen,
                              bmaxlen);

  GT_FREEARRAY(&mlist, GtSeedExtendSeedPair);
  if (benchmark)
    gt_timer_delete(timer);
}
