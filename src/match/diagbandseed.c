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
#include "core/codetype.h"
#include "core/complement.h"
#include "core/encseq.h"
#include "core/minmax.h"
#include "core/radix_sort.h"
#include "core/timer_api.h"
#include "core/unused_api.h"
#include "match/diagbandseed.h"
#include "match/ft-front-prune.h"
#include "match/kmercodes.h"
#include "match/seed-extend.h"
#include "match/sfx-mappedstr.h"
#include "match/sfx-suffixer.h"

#define GT_DIAGBANDSEED_SEQNUM_UNDEF UINT_MAX
#define DEBUG_SEEDPAIR
#undef  DEBUG_SEED_REPORT
#undef  DEBUG_GET_KMERS

typedef uint32_t GtDiagbandseedPosition;
typedef uint32_t GtDiagbandseedSeqnum;
typedef uint32_t GtDiagbandseedScore;
typedef struct GtDiagbandseedProcKmerInfo GtDiagbandseedProcKmerInfo;

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

struct GtDiagbandseedProcKmerInfo {
  GtDiagbandseedKmerPos *list;
  GtUword numberofkmerscollected;
  GtDiagbandseedSeqnum seqnum;
  GtDiagbandseedPosition endpos;
  const GtEncseq *encseq;
  unsigned int kmerlen;
  GtReadmode readmode;
  bool has_short_sequences;
};

/* Add given code and its seqnum and position to a list. */
static void gt_diagbandseed_processkmercode(void *prockmerinfo,
                                            bool firstinrange,
                                            GtUword position,
                                            GtCodetype code)
{
  GtDiagbandseedProcKmerInfo *arg = (GtDiagbandseedProcKmerInfo *) prockmerinfo;
  GtDiagbandseedKmerPos *kmerposptr = arg->list + arg->numberofkmerscollected;
  if (firstinrange == true) {
    arg->endpos = (GtDiagbandseedPosition)(arg->kmerlen - 1);
    if (arg->has_short_sequences) {
      arg->seqnum = (GtDiagbandseedSeqnum)gt_encseq_seqnum(arg->encseq,
                                                           position);
    } else if (arg->seqnum == GT_DIAGBANDSEED_SEQNUM_UNDEF) {
      arg->seqnum = 0;
    } else {
      arg->seqnum++;
    }
  }
  if (arg->readmode == GT_READMODE_FORWARD) {
    kmerposptr->code = code;
  } else {
    kmerposptr->code = gt_kmercode_reverse(code, arg->kmerlen);
  }
  kmerposptr->seqnum = arg->seqnum;
  kmerposptr->endpos = arg->endpos++;
  arg->numberofkmerscollected++;
}

/* Uses GtKmercodeiterator for fetching the kmers. */
static void gt_diagbandseed_get_kmers_kciter(GtDiagbandseedProcKmerInfo *pkinfo)
{
  GtKmercodeiterator *kc_iter = NULL;
  const GtKmercode *kmercode = NULL;
  bool firstinrange = true;
  GtDiagbandseedPosition position;
  gt_assert(pkinfo != NULL);

  kc_iter = gt_kmercodeiterator_encseq_new(pkinfo->encseq, pkinfo->readmode,
                                           pkinfo->kmerlen, 0);
  while ((kmercode = gt_kmercodeiterator_encseq_next(kc_iter)) != NULL) {
    if (!kmercode->definedspecialposition) {
      position = gt_kmercodeiterator_encseq_get_currentpos(kc_iter) - 1;
      gt_diagbandseed_processkmercode((void *)pkinfo, firstinrange, position,
                                      kmercode->code);
      firstinrange = false;
    } else {
      firstinrange = true;
      /* maybe allow comparison of N-containing kmers?! */
    }
  }
  gt_kmercodeiterator_delete(kc_iter);
}

/* Returns a GtDiagbandseedKmerPos list of k-mers from a given encseq. */
GtUword gt_diagbandseed_get_kmers(GtDiagbandseedKmerPos *list,
                                  const GtEncseq *encseq,
                                  unsigned int kmerlen,
                                  GtReadmode readmode)
{
  GtDiagbandseedProcKmerInfo pkinfo;
  gt_assert(encseq != NULL && list != NULL);

  pkinfo.list = list;
  pkinfo.numberofkmerscollected = 0;
  pkinfo.seqnum = pkinfo.endpos = GT_DIAGBANDSEED_SEQNUM_UNDEF;
  pkinfo.encseq = encseq;
  pkinfo.kmerlen = kmerlen;
  pkinfo.readmode = readmode;
  pkinfo.has_short_sequences = (gt_encseq_min_seq_length(encseq) <
                                (GtUword)kmerlen) ? true : false;

  if (gt_encseq_has_twobitencoding(encseq) && gt_encseq_wildcards(encseq) == 0)
  {
    /* Use fast access to encseq, requires 2bit-enc and absence of wildcards. */
    getencseqkmers_twobitencoding(encseq,
                                  readmode,
                                  kmerlen,
                                  kmerlen,
                                  false,
                                  gt_diagbandseed_processkmercode,
                                  (void *) &pkinfo,
                                  NULL,
                                  NULL);
  } else {
    /* Use GtKmercodeiterator for encseq access */
    gt_diagbandseed_get_kmers_kciter(&pkinfo);
  }
  return pkinfo.numberofkmerscollected;
}

/* Returns a GtDiagbandseedSeedPair list of equal kmers from lists a and b. */
void gt_diagbandseed_merge(GtArrayGtDiagbandseedSeedPair *mlist,
                           const GtDiagbandseedKmerPos *alist, GtUword alen,
                           const GtDiagbandseedKmerPos *blist, GtUword blen,
                           unsigned int endposdiff, GtUword maxfreq,
                           bool selfcomp)
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
            if (!selfcomp || aiter->seqnum < biter->seqnum ||
                (aiter->seqnum == biter->seqnum && aiter->endpos + endposdiff <
                 biter->endpos)) {
              /* no duplicates from the same dataset */
              GtDiagbandseedSeedPair *seedptr = NULL;
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

int gt_diagbandseed_process_seeds(const GtEncseq *aencseq,
                                  const GtEncseq *bencseq,
                                  const GtDiagbandseed *arg,
                                  const GtArrayGtDiagbandseedSeedPair *mlist,
                                  GtUword amaxlen,
                                  GtUword bmaxlen,
                                  GtError *err)
{
  int had_err = 0;
  const bool selfcomp = (bencseq == aencseq) ? true : false;
  const GtUword mlen = mlist->nextfreeGtDiagbandseedSeedPair; /* mlist length */
  const GtDiagbandseedSeedPair *lm = mlist->spaceGtDiagbandseedSeedPair;
  const GtUword ndiags = (amaxlen >> arg->logdiagbandwidth) +
                         (bmaxlen >> arg->logdiagbandwidth) + 2;
  /* minimum segment length */
  const GtUword minhit = (arg->mincoverage-1) / arg->seedlength + 1;
  GtUword diag, idx, maxsegm, nextsegm = 0;
  GtDiagbandseedScore *score = NULL;
  GtDiagbandseedPosition *lastp = NULL;

  if (mlen < minhit)
    return 0;
  if (!selfcomp &&
      (arg->extendxdropinfo != NULL || arg->extendgreedyinfo != NULL)) {
    gt_error_set(err, "comparison of two encseqs not implemented");
    return -1;
  }
  score = gt_calloc(ndiags, sizeof *score);
  lastp = gt_calloc(ndiags, sizeof *lastp);
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
             >> arg->logdiagbandwidth;
      if (lm[nextsegm].bpos >= arg->seedlength + lastp[diag]) {
        score[diag] += arg->seedlength;
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
             >> arg->logdiagbandwidth;
      if ((GtUword)(MAX(score[diag+1], score[diag-1])) + (GtUword)(score[diag])
          >= arg->mincoverage)
      {
        GtUword astart = lm[idx].apos + gt_encseq_seqstartpos
          (aencseq, lm[idx].aseqnum) + 1 - arg->seedlength;
        GtUword bstart = lm[idx].bpos + gt_encseq_seqstartpos
          (bencseq, lm[idx].bseqnum) + 1 - arg->seedlength;
#ifdef DEBUG_SEED_REPORT
        printf("SP(%d,%d,%d,%d) ",
               lm[idx].aseqnum, lm[idx].bseqnum, lm[idx].apos, lm[idx].bpos);
#endif
        if (arg->extendgreedyinfo != NULL) {
          had_err = gt_simplegreedyselfmatchoutput((void *)
                                                   arg->extendgreedyinfo,
                                                   aencseq, arg->seedlength,
                                                   astart, bstart, err);
        } else if (arg->extendxdropinfo != NULL) {
          had_err = gt_simplexdropselfmatchoutput((void *)arg->extendxdropinfo,
                                                  aencseq, arg->seedlength,
                                                  astart, bstart, err);
        }
        /* else: no seed extension */
      }
    }

    /* reset diagonal band scores */
    for (idx = currsegm; idx < nextsegm; idx++) {
      diag = (amaxlen + (GtUword)lm[idx].bpos - (GtUword)lm[idx].apos)
             >> arg->logdiagbandwidth;
      score[diag] = 0;
      lastp[diag] = 0;
    }
  }
  gt_free(score);
  gt_free(lastp);
  return had_err;
}

int gt_diagbandseed_run(const GtEncseq *aencseq, const GtEncseq *bencseq,
                        const GtDiagbandseed *arg, GtError *err)
{
  GtDiagbandseedKmerPos *alist = NULL, *blist = NULL;
  GtArrayGtDiagbandseedSeedPair mlist;
  GtRadixsortinfo* rdxinfo = NULL;
  GtUword alen, blen, mlen;
  const unsigned int kmerlen = arg->seedlength;
  const unsigned int endposdiff = arg->overlappingseeds ? 0 : kmerlen - 1;
  const bool selfcomp = (bencseq == aencseq) ? true : false;
  const GtUword amaxlen = gt_encseq_max_seq_length(aencseq);
  const GtUword bmaxlen = gt_encseq_max_seq_length(bencseq);
  /* estimate number of kmers for alist and blist */
  const GtUword ankmers = gt_encseq_total_length(aencseq) -
                          MIN(kmerlen-1, gt_encseq_min_seq_length(aencseq)) *
                          gt_encseq_num_of_sequences(aencseq);
  const GtUword bnkmers = gt_encseq_total_length(bencseq) -
                          MIN(kmerlen-1, gt_encseq_min_seq_length(bencseq)) *
                          gt_encseq_num_of_sequences(bencseq);
  GtTimer *btimer = NULL, *vtimer = NULL;
  int had_err = 0;

  if (amaxlen < kmerlen || bmaxlen < kmerlen) {
    gt_error_set(err,
                 "maximum sequence length too short for required seedlength");
    if (arg->benchmark) {
      printf("0.000000\n");
    }
    return -1;
  }

  if (arg->benchmark) {
    btimer = gt_timer_new();
    gt_timer_start(btimer);
  }

  if (arg->verbose) {
    vtimer = gt_timer_new();
    printf("Start fetching (at most "GT_WU") k-mers for list A...\n", ankmers);
    gt_timer_start(vtimer);
  }

  /* prepare list of kmers from aencseq and sort */
  alist = gt_malloc(ankmers * sizeof *alist);
  alen = gt_diagbandseed_get_kmers(alist, aencseq, kmerlen,
                                   GT_READMODE_FORWARD);

  rdxinfo = gt_radixsort_new_ulongpair(alen);
  gt_radixsort_inplace_GtUwordPair((GtUwordPair*)alist, alen);
  gt_radixsort_delete(rdxinfo);
  if (arg->verbose) {
    printf("Found and sorted "GT_WU" k-mers ", alen);
    gt_timer_show_formatted(vtimer, "in "GT_WD".%06ld seconds.\n", stdout);
  }

  if (!selfcomp || arg->mirror) {
    /* allocate second list */
    if (arg->verbose) {
      printf("Start fetching (at most "GT_WU") k-mers for list B...\n",
             arg->mirror ? 2 * bnkmers : bnkmers);
      gt_timer_start(vtimer);
    }
    if (arg->mirror) {
      blist = gt_malloc(2 * bnkmers * sizeof *blist);
    } else {
      blist = gt_malloc(bnkmers * sizeof *blist);
    }
    /* fill list with forward kmers */
    if (selfcomp) {
      memcpy(blist, alist, alen * sizeof *alist);
      blen = alen;
    } else {
      blen = gt_diagbandseed_get_kmers(blist, bencseq, kmerlen,
                                       GT_READMODE_FORWARD);
    }
    /* add reverse complement kmers */
    if (arg->mirror) {
      blen += gt_diagbandseed_get_kmers(blist + blen, bencseq, kmerlen,
                                        GT_READMODE_COMPL);
    }

    /* sort blist by kmercode */
    rdxinfo = gt_radixsort_new_ulongpair(blen);
    gt_radixsort_inplace_GtUwordPair((GtUwordPair*)blist, blen);
    gt_radixsort_delete(rdxinfo);
    if (arg->verbose) {
      printf("Found and sorted "GT_WU" k-mers ", blen);
      gt_timer_show_formatted(vtimer, "in "GT_WD".%06ld seconds.\n", stdout);
    }
  } else {
    /* compare reads of encseq A with themselves */
    blist = alist;
    blen = alen;
  }

#ifdef DEBUG_GET_KMERS
  for (GtDiagbandseedKmerPos *a = alist; a < alist+alen; a++) {
    printf("a) Kmer (%lX,%d,%d)\n", a->code, a->endpos, a->seqnum);
  }
  for (GtDiagbandseedKmerPos *b = blist; b < blist+blen; b++) {
    printf("b) Kmer (%lX,%d,%d)\n", b->code, b->endpos, b->seqnum);
  }
#endif

  if (arg->verbose) {
    printf("Start building seed pairs on equal k-mers...\n");
    gt_timer_start(vtimer);
  }
  /* create mlist of SeedPairs */
  GT_INITARRAY(&mlist,GtDiagbandseedSeedPair);
  gt_diagbandseed_merge(&mlist, alist, alen, blist, blen, endposdiff,
                        arg->maxfreq, selfcomp);
  gt_free(alist);
  if (!selfcomp || arg->mirror)
    gt_free(blist);

  /* sort mlist */
  mlen = mlist.nextfreeGtDiagbandseedSeedPair;
  rdxinfo = gt_radixsort_new_uint64keypair(mlen);
  gt_radixsort_inplace_Gtuint64keyPair((Gtuint64keyPair*)mlist.
                                       spaceGtDiagbandseedSeedPair, mlen);
  gt_radixsort_delete(rdxinfo);
  if (arg->verbose) {
    printf("Collected and sorted "GT_WU" seed pairs ", mlen);
    gt_timer_show_formatted(vtimer, "in "GT_WD".%06ld seconds.\n", stdout);
  }

  /* verify SeedPairs in the sequences */
  if (arg->verify && mlen != 0) {
    GtDiagbandseedSeedPair *j = mlist.spaceGtDiagbandseedSeedPair;
    GtDiagbandseedSeedPair *last = j + mlen;
    char *buf1 = gt_malloc(3 * (1 + kmerlen) * sizeof *buf1);
    char *buf2 = buf1 + 1 + kmerlen;
    char *buf3 = buf2 + 1 + kmerlen;
    while (j < last) {
      char *idx = NULL;
      GtDiagbandseedPosition a = j->apos + gt_encseq_seqstartpos(aencseq,
                                                                 j->aseqnum);
      GtDiagbandseedPosition b = j->bpos + gt_encseq_seqstartpos(bencseq,
                                                                 j->bseqnum);
      gt_encseq_extract_decoded(aencseq, buf1, a + 1 - kmerlen, a);
      gt_encseq_extract_decoded(bencseq, buf2, b + 1 - kmerlen, b);
      buf1[kmerlen] = buf2[kmerlen] = buf3[kmerlen] = '\0';
#ifdef DEBUG_SEEDPAIR
      printf("SeedPair (%d,%d,%d,%d)\n", j->aseqnum, j->bseqnum, j->apos,
             j->bpos);
#endif
      for (idx = buf3; idx < buf3 + kmerlen; idx++)
        gt_complement(idx, buf2[kmerlen + buf3 - idx - 1], NULL);
      if (strcmp(buf1, buf2) != 0 && (!arg->mirror || strcmp(buf1, buf3) != 0))
      {
        gt_error_set(err, "wrong seed(%d,%d,%d,%d): %s != %s / %s\n",
                     j->aseqnum, j->bseqnum, j->apos, j->bpos, buf1, buf2,
                     buf3);
        had_err = -1;
      }
      j++;
    }
    gt_free(buf1);
  }

  /* process SeedPairs */
  if (had_err == 0 && arg->verbose) {
    printf("Start seed pair extension...\n");
    gt_timer_start(vtimer);
  }
  if (had_err == 0 && mlen != 0)
    had_err = gt_diagbandseed_process_seeds(aencseq, bencseq, arg, &mlist,
                                            amaxlen, bmaxlen, err);
  GT_FREEARRAY(&mlist, GtDiagbandseedSeedPair);

  if (arg->verbose) {
    if (had_err == 0)
      gt_timer_show_formatted(vtimer, "Completed extension in "GT_WD".%06ld "
                              "seconds.\n", stdout);
    gt_timer_stop(vtimer);
    gt_timer_delete(vtimer);
  }

  if (arg->benchmark) {
    gt_timer_show_formatted(btimer, "Overall time "GT_WD".%06ld seconds.\n",
                            stdout);
    gt_timer_delete(btimer);
  }
  return had_err;
}
