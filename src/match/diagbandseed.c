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
#include "match/diagbandseed.h"
#include "match/sfx-mappedstr.h"

#define DEBUG_SEEDPAIR
#undef  DEBUG_SEED_REPORT

typedef uint32_t GtDiagbandseedPosition;
typedef uint32_t GtDiagbandseedSeqnum;
typedef uint32_t GtDiagbandseedScore;

struct GtDiagbandseedKmerPos {
  GtCodetype code;            /* only sort criterion */
  GtDiagbandseedSeqnum seqnum;
  GtDiagbandseedPosition endpos;
};

struct GtDiagbandseedSeedPair {
  GtDiagbandseedSeqnum bseqnum; /*  2nd important sort criterion */
  GtDiagbandseedSeqnum aseqnum; /* most important sort criterion */
  GtDiagbandseedPosition apos;
  GtDiagbandseedPosition bpos;  /*  3rd important sort criterion */
};

/* Returns a GtDiagbandseedKmerPos list of k-mers from a given encseq. */
GtUword gt_diagbandseed_get_kmers(GtDiagbandseedKmerPos *list,
                                  const GtEncseq *encseq, unsigned int kmerlen)
{
  const GtKmercode *kmercode;
  GtKmercodeiterator *kc_iter;
  GtUword numberofkmerscollected = 0;
  bool nextsegm = true;
  GtUword seqnum = 0, endpos = 0;

  gt_assert(list != NULL && encseq != NULL);
  kc_iter = gt_kmercodeiterator_encseq_new(encseq, GT_READMODE_FORWARD, kmerlen,
                                           0);
  while ((kmercode = gt_kmercodeiterator_encseq_next(kc_iter)) != NULL) {
    if (!kmercode->definedspecialposition) {
      /* store (kmercode, seqnum, endpos) in array */
      if (nextsegm == true) {
        endpos = gt_kmercodeiterator_encseq_get_currentpos(kc_iter) - 1;
        seqnum = gt_encseq_seqnum(encseq, endpos);
        endpos = kmerlen - 1;
        nextsegm = false;
      } else {
        endpos++;
      }
      GtDiagbandseedKmerPos *kmerposptr = list + numberofkmerscollected;
      kmerposptr->code = kmercode->code;
      kmerposptr->endpos = endpos;
      kmerposptr->seqnum = seqnum;
      numberofkmerscollected++;
    } else {
      nextsegm = true;
      /* TODO: allow comparison of N-containing kmers */
    }
  }
  gt_kmercodeiterator_delete(kc_iter);
  return numberofkmerscollected;
}

/* Returns a GtDiagbandseedSeedPair list of equal kmers from lists a and b. */
void gt_diagbandseed_merge(GtArrayGtDiagbandseedSeedPair *mlist,
                           const GtDiagbandseedKmerPos *alist, GtUword alen,
                           const GtDiagbandseedKmerPos *blist, GtUword blen,
                           unsigned int kmerlen, GtUword maxfreq)
{
  const GtDiagbandseedKmerPos *aptr = alist, *bptr = blist, *aend, *bend;
  const GtUword array_incr = 256;

  gt_assert(alist != NULL && blist != NULL && mlist != NULL);
  aend = aptr + alen;
  bend = bptr + blen;
  while (aptr < aend && bptr < bend) {
    if (aptr->code < bptr->code) {
      aptr++;
    } else if (aptr->code > bptr->code) {
      bptr++;
    } else {
      /* equality: count frequency of current k-mer in both lists */
      const GtDiagbandseedKmerPos *aiter, *biter;
      for (aiter = aptr; aiter < aend && aiter->code == bptr->code; aiter++) {
        /* nothing */
      }
      for (biter = bptr; biter < bend && biter->code == aptr->code; biter++) {
        /* nothing */
      }
      if ((GtUword)(aiter - aptr) <= maxfreq &&
          (GtUword)(biter - bptr) <= maxfreq) {
        /* add all equal k-mers */
        const GtDiagbandseedKmerPos *asegm_end = aiter, *bsegm_end = biter;
        for (aiter = aptr; aiter < asegm_end; aiter++) {
          for (biter = bptr; biter < bsegm_end; biter++) {
            if (alist != blist || aiter->seqnum < biter->seqnum ||
                (aiter->seqnum == biter->seqnum && aiter->endpos + kmerlen <=
                 biter->endpos)) {
              /* no duplicates from the same dataset */
              GtDiagbandseedSeedPair *seedptr;
              GT_GETNEXTFREEINARRAY(seedptr, mlist, GtDiagbandseedSeedPair,
                                    array_incr + 0.2 *
                                    mlist->allocatedGtDiagbandseedSeedPair);
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

bool gt_diagbandseed_is_seed(GT_UNUSED const GtDiagbandseedSeedPair *entry,
                             const GtDiagbandseedScore *score,
                             GtUword mincoverage, GtUword diag)
{
  /* TODO: add path criterion */
  gt_assert(score != NULL);
  if ((GtUword)(MAX(score[diag+1], score[diag-1])) + (GtUword)(score[diag])
      >= mincoverage) {
    return true;
  } else {
    return false;
  }
}

void gt_diagbandseed_find_seeds(const GtArrayGtDiagbandseedSeedPair *mlist,
                                unsigned int kmerlen, GtUword mincoverage,
                                GtUword log_diagbandwidth, GtUword amaxlen,
                                GtUword bmaxlen)
{
  const GtUword mlen = mlist->nextfreeGtDiagbandseedSeedPair; /* mlist length */
  const GtDiagbandseedSeedPair *lm = mlist->spaceGtDiagbandseedSeedPair;
  const GtUword ndiags = (amaxlen >> log_diagbandwidth) +
                         (bmaxlen >> log_diagbandwidth) + 2;
  const GtUword minhit = (mincoverage-1) / kmerlen + 1; /* min segment length */
  GtUword diag, idx, maxsegm, nextsegm = 0;
  GtDiagbandseedScore *score = gt_calloc(ndiags, sizeof *score);
  GtDiagbandseedPosition *lastp = gt_calloc(ndiags, sizeof *lastp);

  if (mlen < minhit)
    return;
  maxsegm = mlen - minhit;

  while (nextsegm <= maxsegm) {
    const GtUword currsegm = nextsegm;
    const GtDiagbandseedSeqnum currsegm_aseqnum = lm[currsegm].aseqnum;
    const GtDiagbandseedSeqnum currsegm_bseqnum = lm[currsegm].bseqnum;

    /* if insuffienct number of kmers in segment: skip whole segment */
    if (currsegm_aseqnum != lm[currsegm + minhit - 1].aseqnum ||
        currsegm_bseqnum != lm[currsegm + minhit - 1].bseqnum) {
      do {
        nextsegm++;
      } while (nextsegm < mlen && lm[nextsegm].aseqnum == currsegm_aseqnum &&
               lm[nextsegm].bseqnum == currsegm_bseqnum);
      continue;
    }

    /* calculate diagonal band scores */
    do {
      gt_assert(lm[nextsegm].bpos <= bmaxlen && lm[nextsegm].apos <= amaxlen);
      diag = (amaxlen + (GtUword)lm[nextsegm].bpos - (GtUword)lm[nextsegm].apos)
             >> log_diagbandwidth;
      if (lm[nextsegm].bpos >= kmerlen + lastp[diag]) {
        score[diag] += kmerlen;
      } else {
        gt_assert(lastp[diag] <= lm[nextsegm].bpos);/*if fail: sorted by bpos?*/
        score[diag] = score[diag] + lm[nextsegm].bpos - lastp[diag];
      }
      lastp[diag] = lm[nextsegm].bpos;
      nextsegm++;
    } while (nextsegm < mlen && lm[nextsegm].aseqnum == currsegm_aseqnum &&
             lm[nextsegm].bseqnum == currsegm_bseqnum);

    /* report seeds */
    for (idx = currsegm; idx < nextsegm; idx++) {
      gt_assert(lm[idx].apos <= amaxlen);
      diag = (amaxlen + (GtUword)lm[idx].bpos - (GtUword)lm[idx].apos)
             >> log_diagbandwidth;
      if (gt_diagbandseed_is_seed(&lm[idx], score, mincoverage, diag)) {
#ifdef DEBUG_SEED_REPORT
        printf("report SeedPair (%d,%d,%d,%d), score[%lu]=%d\n",
               lm[idx].aseqnum, lm[idx].bseqnum, lm[idx].apos, lm[idx].bpos,
               diag, MAX(score[diag+1], score[diag-1]) + score[diag]);
#endif
      }
    }

    /* reset diagonal band scores */
    for (idx = currsegm; idx < nextsegm; idx++) {
      diag = (amaxlen + (GtUword)lm[idx].bpos - (GtUword)lm[idx].apos)
             >> log_diagbandwidth;
      score[diag] = 0;
      lastp[diag] = 0;
    }
  }
  gt_free(score);
  gt_free(lastp);
}

void gt_diagbandseed_run(const GtEncseq *aencseq, const GtEncseq *bencseq,
                         const GtDiagbandseed *arg)
{
  GtDiagbandseedKmerPos *alist, *blist;
  GtArrayGtDiagbandseedSeedPair mlist;
  GtRadixsortinfo* rdxinfo;
  GtUword alen, blen;
  const unsigned int kmerlen = arg->dbs_seedlength;
  const bool two_files = (bencseq != aencseq) ? true : false;
  const GtUword amaxlen = gt_encseq_max_seq_length(aencseq);
  const GtUword bmaxlen = gt_encseq_max_seq_length(bencseq);
  const GtUword ankmers = gt_encseq_total_length(aencseq) -
                          MIN(kmerlen-1, gt_encseq_min_seq_length(aencseq)) *
                          gt_encseq_num_of_sequences(aencseq);
  const GtUword bnkmers = gt_encseq_total_length(bencseq) -
                          MIN(kmerlen-1, gt_encseq_min_seq_length(bencseq)) *
                          gt_encseq_num_of_sequences(bencseq);
  GtTimer *timer;

  if (amaxlen < kmerlen || bmaxlen < kmerlen) {
    /*printf("maximum sequence length too short\n");*/
    if (arg->benchmark) {
      printf("0.000000,0,0\n");
    }
    return;
  }

  alist = gt_malloc(ankmers * sizeof *alist);
  alen = gt_diagbandseed_get_kmers(alist, aencseq, kmerlen);
  if (arg->benchmark) {
    timer = gt_timer_new();
    gt_timer_start(timer);
  }
  rdxinfo = gt_radixsort_new_ulongpair(alen);
  gt_radixsort_inplace_GtUwordPair((GtUwordPair*)alist, alen);
  gt_radixsort_delete(rdxinfo);

  if (two_files) {
    blist = gt_malloc(bnkmers * sizeof *blist);
    blen = gt_diagbandseed_get_kmers(blist, bencseq, kmerlen);
    rdxinfo = gt_radixsort_new_ulongpair(blen);
    gt_radixsort_inplace_GtUwordPair((GtUwordPair*)blist, blen);
    gt_radixsort_delete(rdxinfo);
  } else {
    /* compare all reads of encseq A with themselves */
    blist = alist;
    blen = alen;
  }

  GT_INITARRAY(&mlist,GtDiagbandseedSeedPair);
  gt_diagbandseed_merge(&mlist, alist, alen, blist, blen, kmerlen,
                        arg->dbs_maxfreq);
  gt_free(alist);
  if (two_files)
    gt_free(blist);
  rdxinfo = gt_radixsort_new_uint64keypair(mlist.
                                           nextfreeGtDiagbandseedSeedPair);
  gt_radixsort_inplace_Gtuint64keyPair((Gtuint64keyPair*)mlist.
                                       spaceGtDiagbandseedSeedPair,
                                       mlist.nextfreeGtDiagbandseedSeedPair);
  gt_radixsort_delete(rdxinfo);
  if (arg->benchmark) {
    gt_timer_stop(timer);
    gt_timer_show_formatted(timer, GT_WD ".%06ld,"GT_WD","GT_WD"\n", stdout);
  }

  if (arg->verify && mlist.nextfreeGtDiagbandseedSeedPair != 0) {
    GtDiagbandseedSeedPair *j = mlist.spaceGtDiagbandseedSeedPair;
    GtDiagbandseedSeedPair *last = j + mlist.nextfreeGtDiagbandseedSeedPair;
    char *buf1 = gt_malloc(1 + kmerlen * sizeof *buf1);
    char *buf2 = gt_malloc(1 + kmerlen * sizeof *buf2);
    while (j < last) {
      GtDiagbandseedPosition a = j->apos + gt_encseq_seqstartpos(aencseq,
                                                                 j->aseqnum);
      GtDiagbandseedPosition b = j->bpos + gt_encseq_seqstartpos(bencseq,
                                                                 j->bseqnum);
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
        gt_assert(false);
      }
      j++;
    }
    gt_free(buf1);
    gt_free(buf2);
  }

  if (mlist.nextfreeGtDiagbandseedSeedPair != 0)
    gt_diagbandseed_find_seeds(&mlist, kmerlen, arg->dbs_mincoverage,
                               arg->dbs_logdiagbandwidth, amaxlen, bmaxlen);

  GT_FREEARRAY(&mlist, GtDiagbandseedSeedPair);
  if (arg->benchmark)
    gt_timer_delete(timer);
}
