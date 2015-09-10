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

#include <limits.h>
#include <stdio.h>
#include <string.h>
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
#include "match/querymatch.h"
#include "match/querymatch-align.h"
#include "match/seed-extend.h"
#include "match/sfx-mappedstr.h"
#include "match/sfx-suffixer.h"

#define GT_DIAGBANDSEED_SEQNUM_UNDEF UINT_MAX

typedef uint32_t GtDiagbandseedPosition;
typedef uint32_t GtDiagbandseedSeqnum;
typedef uint32_t GtDiagbandseedScore;
typedef struct GtDiagbandseedProcKmerInfo GtDiagbandseedProcKmerInfo;
typedef const GtQuerymatch *(*GtDiagbandseedExtendFunc)(void *,
                                                        const GtEncseq *,
                                                        GtUword,
                                                        GtUword,
                                                        GtUword,
                                                        GtUword,
                                                        GtUword);

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
  unsigned int seedlength;
  GtReadmode readmode;
  bool has_short_sequences;
};

/* Add given code and its seqnum and position to a list. */
static void gt_diagbandseed_processkmercode(void *prockmerinfo,
                                            bool firstinrange,
                                            GtUword position,
                                            GtCodetype code)
{
  GtDiagbandseedProcKmerInfo *arg;
  GtDiagbandseedKmerPos *kmerposptr;
  GtUword seqstartpos = 0;

  gt_assert(prockmerinfo != NULL);
  arg = (GtDiagbandseedProcKmerInfo *) prockmerinfo;
  kmerposptr = arg->list + arg->numberofkmerscollected;

  /* save seqnum and, if necessary, reset endpos */
  if (firstinrange == true) {
    if (arg->has_short_sequences) {
      arg->seqnum = (GtDiagbandseedSeqnum)gt_encseq_seqnum(arg->encseq,
                                                           position);
    } else if (arg->seqnum == GT_DIAGBANDSEED_SEQNUM_UNDEF) {
      arg->seqnum = 0;
    } else {
      gt_assert(arg->seqnum != UINT_MAX);
      arg->seqnum++;
    }
    seqstartpos = gt_encseq_seqstartpos(arg->encseq, arg->seqnum);
    gt_assert(position >= seqstartpos);
    arg->endpos = (GtDiagbandseedPosition)(position - seqstartpos);
  }
  kmerposptr->seqnum = arg->seqnum;

  /* save k-mer code */
  if (arg->readmode == GT_READMODE_FORWARD) {
    kmerposptr->code = code;
  } else {
    kmerposptr->code = gt_kmercode_reverse(code, arg->seedlength);
  }

  /* save endpos */
  gt_assert(arg->endpos != UINT_MAX);
  kmerposptr->endpos = arg->endpos++;
  arg->numberofkmerscollected++;
}

/* Uses GtKmercodeiterator for fetching the kmers. */
static void gt_diagbandseed_get_kmers_kciter(GtDiagbandseedProcKmerInfo *pkinfo)
{
  GtKmercodeiterator *kc_iter = NULL;
  const GtKmercode *kmercode = NULL;
  GtDiagbandseedPosition position = 0;
  bool firstinrange = true;

  /* initialise GtKmercodeiterator */
  gt_assert(pkinfo != NULL);
  kc_iter = gt_kmercodeiterator_encseq_new(pkinfo->encseq,
                                           pkinfo->readmode,
                                           pkinfo->seedlength,
                                           0);

  while ((kmercode = gt_kmercodeiterator_encseq_next(kc_iter)) != NULL) {
    if (!kmercode->definedspecialposition) {
      position = gt_kmercodeiterator_encseq_get_currentpos(kc_iter) - 1;
      gt_diagbandseed_processkmercode((void *)pkinfo,
                                      firstinrange,
                                      position,
                                      kmercode->code);
      firstinrange = false;
    } else {
      /* separator: next sequence */
      firstinrange = true;
    }
  }
  gt_kmercodeiterator_delete(kc_iter);
}

/* Returns a GtDiagbandseedKmerPos list of k-mers from a given encseq. */
GtUword gt_diagbandseed_get_kmers(GtDiagbandseedKmerPos *list,
                                  const GtEncseq *encseq,
                                  unsigned int seedlength,
                                  GtReadmode readmode)
{
  GtDiagbandseedProcKmerInfo pkinfo;

  gt_assert(list != NULL);
  pkinfo.list = list;
  pkinfo.numberofkmerscollected = 0;
  pkinfo.seqnum = pkinfo.endpos = GT_DIAGBANDSEED_SEQNUM_UNDEF;
  gt_assert(encseq != NULL);
  pkinfo.encseq = encseq;
  pkinfo.seedlength = seedlength;
  pkinfo.readmode = readmode;
  pkinfo.has_short_sequences = (gt_encseq_wildcards(encseq) != 0 ||
                                gt_encseq_min_seq_length(encseq) <
                                (GtUword)seedlength) ? true : false;

  if (gt_encseq_has_twobitencoding(encseq) && gt_encseq_wildcards(encseq) == 0)
  {
    /* Use fast access to encseq, requires 2bit-enc and absence of wildcards. */
    getencseqkmers_twobitencoding(encseq,
                                  readmode,
                                  seedlength,
                                  seedlength,
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
                           GtUword *maxfreq,
                           GtUword maxgram,
                           GtUword memlimit,
                           GtUword *histogram,
                           unsigned int endposdiff,
                           bool selfcomp)
{
  const GtDiagbandseedKmerPos *aptr = alist, *bptr = blist, *aend, *bend;
  const GtUword array_incr = 256;
  GtUword frequency = 0;

  gt_assert(alist != NULL && blist != NULL && mlist != NULL && maxfreq != NULL);
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
      frequency = (GtUword) MAX(aiter - aptr, biter - bptr);
      if (frequency <= *maxfreq) {
        /* add all equal k-mers */
        const GtDiagbandseedKmerPos *asegm_end = aiter, *bsegm_end = biter;
        frequency = MIN(maxgram, frequency);
        gt_assert(frequency > 0);
        for (aiter = aptr; aiter < asegm_end; aiter++) {
          for (biter = bptr; biter < bsegm_end; biter++) {
            if (!selfcomp ||
                aiter->seqnum < biter->seqnum ||
                (aiter->seqnum == biter->seqnum &&
                 aiter->endpos + endposdiff <= biter->endpos)) {
              /* no duplicates from the same dataset */
              if (histogram == NULL) {
                /* save SeedPair in mlist */
                GtDiagbandseedSeedPair *seedptr = NULL;
                GT_GETNEXTFREEINARRAY(seedptr,
                                      mlist,
                                      GtDiagbandseedSeedPair,
                                      array_incr + 0.2 *
                                      mlist->allocatedGtDiagbandseedSeedPair);
                seedptr->bseqnum = biter->seqnum;
                seedptr->aseqnum = aiter->seqnum;
                seedptr->bpos = biter->endpos;
                seedptr->apos = aiter->endpos;
              } else {
                /* count seed pair frequency in histogram */
                histogram[frequency - 1]++;
              }
            }
          }
        }
      } /* else: ignore all equal elements */
      aptr = aiter;
      bptr = biter;
    }
  }
  if (histogram != NULL) {
    /* calculate available memory, take 98% of memlimit */
    GtUword count = 0, mem_used = 0, mem_avail = 0.98 * memlimit;
    if (alist == blist) {
      gt_assert(selfcomp);
      mem_used = alen * sizeof *alist;
    } else {
      mem_used = (alen + blen) * sizeof *alist;
    }
    if (mem_avail > mem_used) {
      mem_avail = (mem_avail - mem_used) / sizeof (GtDiagbandseedSeedPair);
    } else {
      mem_avail = 0;
      *maxfreq = 0;
    }

    /* there is enough free memory */
    if (mem_avail > 0) {
      /* count seed pairs until available memory reached */
      for (frequency = 1; frequency <= maxgram && count < mem_avail;
           frequency++) {
        count += histogram[frequency - 1];
      }
      if (count > mem_avail) {
        gt_assert(frequency >= 2 && count >= histogram[frequency]);
        frequency -= 2;
        count -= histogram[frequency];
      } else if (frequency == maxgram + 1) {
        frequency = GT_UWORD_MAX;
      }
      *maxfreq = MIN(*maxfreq, frequency);
    }

    /* determine minimum required memory for error message */
    if (*maxfreq <= 1 && selfcomp) {
      count = (histogram[0] + histogram[1]) * sizeof (GtDiagbandseedSeedPair);
      count = (count + mem_used) / 0.98;
    } else if (*maxfreq == 0) {
      count = histogram[0] * sizeof (GtDiagbandseedSeedPair);
      count = (count + mem_used) / 0.98;
    }
    histogram[maxgram] = count;
  }
}

/* start seed extension for seed pairs in mlist */
int gt_diagbandseed_process_seeds(const GtEncseq *aencseq,
                                  const GtEncseq *bencseq,
                                  const GtArrayGtDiagbandseedSeedPair *mlist,
                                  GtGreedyextendmatchinfo *extendgreedyinfo,
                                  GtXdropmatchinfo *extendxdropinfo,
                                  GtQuerymatchoutoptions *querymatchoutopt,
                                  unsigned int seedlength,
                                  GtUword logdiagbandwidth,
                                  GtUword mincoverage,
                                  GtUword amaxlen,
                                  GtUword bmaxlen,
                                  GtError *err)
{
  GtDiagbandseedScore *score = NULL;
  GtDiagbandseedPosition *lastp = NULL;
  GtDiagbandseedExtendFunc extend_selfmatch_relpos_function = NULL;
  GtProcessinfo_and_querymatchspaceptr info_querymatch;
  const GtDiagbandseedSeedPair *lm = NULL;
  const GtUword ndiags = (amaxlen >> logdiagbandwidth) +
                         (bmaxlen >> logdiagbandwidth) + 2;
  const GtUword minsegmentlen = (mincoverage - 1) / seedlength + 1;
  GtUword mlen = 0, diag = 0, idx = 0, maxsegm = 0, nextsegm = 0;
  int had_err = 0;
  bool firstinrange = true;

  gt_assert(mlist != NULL);
  mlen = mlist->nextfreeGtDiagbandseedSeedPair; /* mlist length  */
  lm = mlist->spaceGtDiagbandseedSeedPair;      /* mlist pointer */

  /* select extension method */
  if (extendgreedyinfo != NULL) {
    info_querymatch.processinfo = (void *)(extendgreedyinfo);
    extend_selfmatch_relpos_function = gt_greedy_extend_selfmatch_relpos;
  } else if (extendxdropinfo != NULL) {
    info_querymatch.processinfo = (void *)(extendxdropinfo);
    extend_selfmatch_relpos_function = gt_xdrop_extend_selfmatch_relpos;
  } else { /* no seed extension */
    return 0;
  }

  if (mlen < minsegmentlen)
    return 0;
  gt_assert(aencseq != NULL && bencseq != NULL);
  if (aencseq != bencseq) {
    gt_error_set(err, "comparison of two encseqs not implemented");
    return -1;
  }

  info_querymatch.querymatchspaceptr = gt_querymatch_new(querymatchoutopt);
  /* score[0] and score[ndiags+1] remain zero for boundary */
  score = gt_calloc(ndiags + 2, sizeof *score);
  lastp = gt_calloc(ndiags, sizeof *lastp);
  maxsegm = mlen - minsegmentlen;

  /* iterate through segments of equal k-mers */
  while (nextsegm <= maxsegm) {
    const GtUword currsegm = nextsegm;
    const GtDiagbandseedSeqnum currsegm_aseqnum = lm[currsegm].aseqnum;
    const GtDiagbandseedSeqnum currsegm_bseqnum = lm[currsegm].bseqnum;

    /* if insuffienct number of kmers in segment: skip whole segment */
    if (currsegm_aseqnum != lm[currsegm + minsegmentlen - 1].aseqnum ||
        currsegm_bseqnum != lm[currsegm + minsegmentlen - 1].bseqnum) {
      do {
        nextsegm++;
      } while (nextsegm < mlen &&
               lm[nextsegm].aseqnum == currsegm_aseqnum &&
               lm[nextsegm].bseqnum == currsegm_bseqnum);
      continue;
    }

    /* calculate diagonal band scores */
    do {
      gt_assert(lm[nextsegm].bpos <= bmaxlen && lm[nextsegm].apos <= amaxlen);
      diag = (amaxlen + (GtUword)lm[nextsegm].bpos - (GtUword)lm[nextsegm].apos)
             >> logdiagbandwidth;
      if (lm[nextsegm].bpos >= seedlength + lastp[diag]) {
        /* no overlap: add seedlength */
        score[diag + 1] += seedlength;
      } else {
        /* overlap: add difference below overlap */
        gt_assert(lastp[diag] <= lm[nextsegm].bpos); /* if fail: sort by bpos */
        score[diag + 1] += lm[nextsegm].bpos - lastp[diag];
      }
      lastp[diag] = lm[nextsegm].bpos;
      nextsegm++;
    } while (nextsegm < mlen && lm[nextsegm].aseqnum == currsegm_aseqnum &&
             lm[nextsegm].bseqnum == currsegm_bseqnum);

    /* test for mincoverage and overlap to previous extension */
    firstinrange = true;
    for (idx = currsegm; idx < nextsegm; idx++) {
      gt_assert(lm[idx].apos <= amaxlen);
      diag = (amaxlen + (GtUword)lm[idx].bpos - (GtUword)lm[idx].apos)
             >> logdiagbandwidth;
      if ((GtUword)MAX(score[diag + 2], score[diag]) + (GtUword)score[diag + 1]
          >= mincoverage)
      {
        /* relative seed start positions */
        GtUword astart = lm[idx].apos + 1 - seedlength;
        GtUword bstart = lm[idx].bpos + 1 - seedlength;

        gt_assert(info_querymatch.querymatchspaceptr != NULL);
        if (firstinrange ||
            gt_querymatch_checkoverlap(info_querymatch.querymatchspaceptr,
                                       lm[idx].bseqnum,
                                       bstart))
        {
          /* extend seed */
          const GtQuerymatch *querymatch;
          querymatch = extend_selfmatch_relpos_function(&info_querymatch,
                                                        aencseq,
                                                        lm[idx].aseqnum,
                                                        astart,
                                                        lm[idx].bseqnum,
                                                        bstart,
                                                        seedlength);
          if (querymatch != NULL) {
            firstinrange = false;
            /* show extension results */
            gt_querymatch_prettyprint(querymatch);
          }
        }
      }
    }

    /* reset diagonal band scores */
    for (idx = currsegm; idx < nextsegm; idx++) {
      diag = (amaxlen + (GtUword)lm[idx].bpos - (GtUword)lm[idx].apos)
             >> logdiagbandwidth;
      score[diag + 1] = 0;
      lastp[diag] = 0;
    }
  }
  gt_querymatch_delete(info_querymatch.querymatchspaceptr);
  gt_free(score);
  gt_free(lastp);
  return had_err;
}

int gt_diagbandseed_run(const GtEncseq *aencseq,
                        const GtEncseq *bencseq,
                        const GtDiagbandseed *arg,
                        GtError *err)
{
  GtDiagbandseedKmerPos *alist = NULL, *blist = NULL;
  GtArrayGtDiagbandseedSeedPair mlist;
  GtRadixsortinfo *rdxinfo = NULL;
  GtTimer *vtimer = NULL;
  GtUword alen = 0, blen = 0, mlen = 0, maxfreq = 0;
  GtUword amaxlen = 0, bmaxlen = 0, ankmers = 0, bnkmers = 0;
  unsigned int endposdiff = 0;
  int had_err = 0;
  const bool selfcomp = (bencseq == aencseq) ? true : false;
  const GtUword maxgram = 10000; /* Cap on k-mer count histogram */

  gt_assert(arg != NULL);
  maxfreq = arg->maxfreq;
  endposdiff = arg->overlappingseeds == false ? arg->seedlength : 1;

  gt_assert(aencseq != NULL && bencseq != NULL);
  amaxlen = gt_encseq_max_seq_length(aencseq);
  bmaxlen = gt_encseq_max_seq_length(bencseq);
  /* estimate number of kmers for alist and blist */
  ankmers = gt_encseq_total_length(aencseq) - MIN(arg->seedlength - 1,
    gt_encseq_min_seq_length(aencseq)) * gt_encseq_num_of_sequences(aencseq);
  bnkmers = gt_encseq_total_length(bencseq) - MIN(arg->seedlength - 1,
    gt_encseq_min_seq_length(bencseq)) * gt_encseq_num_of_sequences(bencseq);

  if (amaxlen < arg->seedlength || bmaxlen < arg->seedlength) {
    gt_error_set(err,
                 "maximum sequence length too short for required seedlength");
    return -1;
  }

  if (arg->verbose) {
    vtimer = gt_timer_new();
    printf("# Start fetching (at most "GT_WU") k-mers for list A...\n",
           ankmers);
    gt_timer_start(vtimer);
  }

  /* prepare list of kmers from aencseq and sort */
  alist = gt_malloc(ankmers * sizeof *alist);
  alen = gt_diagbandseed_get_kmers(alist,
                                   aencseq,
                                   arg->seedlength,
                                   GT_READMODE_FORWARD);

  rdxinfo = gt_radixsort_new_ulongpair(alen);
  gt_radixsort_inplace_GtUwordPair((GtUwordPair *)alist, alen);
  gt_radixsort_delete(rdxinfo);
  if (arg->verbose) {
    printf("# Found and sorted "GT_WU" k-mers ", alen);
    gt_timer_show_formatted(vtimer, "in "GT_WD".%06ld seconds.\n", stdout);
  }

  if (arg->debug_kmer) {
    for (GtDiagbandseedKmerPos *a = alist; a < alist + alen; a++) {
      printf("a) Kmer (%lX,%d,%d)\n", a->code, a->endpos, a->seqnum);
    }
  }

  if (!selfcomp || arg->mirror) {
    /* allocate second list */
    if (arg->verbose) {
      printf("# Start fetching (at most "GT_WU") k-mers for list B...\n",
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
      blen = gt_diagbandseed_get_kmers(blist,
                                       bencseq,
                                       arg->seedlength,
                                       GT_READMODE_FORWARD);
    }
    /* add reverse complement kmers */
    if (arg->mirror) {
      blen += gt_diagbandseed_get_kmers(blist + blen,
                                        bencseq,
                                        arg->seedlength,
                                        GT_READMODE_COMPL);
    }

    /* sort blist by kmercode */
    rdxinfo = gt_radixsort_new_ulongpair(blen);
    gt_radixsort_inplace_GtUwordPair((GtUwordPair *)blist, blen);
    gt_radixsort_delete(rdxinfo);
    if (arg->verbose) {
      printf("# Found and sorted "GT_WU" k-mers ", blen);
      gt_timer_show_formatted(vtimer, "in "GT_WD".%06ld seconds.\n", stdout);
    }

    if (arg->debug_kmer) {
      for (GtDiagbandseedKmerPos *b = blist; b < blist + blen; b++) {
        printf("b) Kmer (%lX,%d,%d)\n", b->code, b->endpos, b->seqnum);
      }
    }
  } else {
    /* compare reads of encseq A with themselves */
    blist = alist;
    blen = alen;
  }

  /* calculate maxfreq from memlimit */
  GT_INITARRAY(&mlist, GtDiagbandseedSeedPair);
  if (arg->memlimit < GT_UWORD_MAX) {
    if (arg->verbose) {
      printf("# Start calculating k-mer frequency histogram...\n");
      gt_timer_start(vtimer);
    }
    GtUword count = 0;
    GtUword *histogram = gt_calloc(maxgram + 1, sizeof *histogram);
    gt_diagbandseed_merge(&mlist,
                          alist,
                          alen,
                          blist,
                          blen,
                          &maxfreq,
                          maxgram,
                          arg->memlimit,
                          histogram,
                          endposdiff,
                          selfcomp);
    count = histogram[maxgram];
    gt_free(histogram);
    if (maxfreq > 1 || (maxfreq == 1 && !selfcomp)) {
      /* allocate mlist according to seed pair count */
      GT_CHECKARRAYSPACEMULTI(&mlist, GtDiagbandseedSeedPair, count);
      if (maxfreq < 10) {
        printf("# Warning: Only k-mers occurring <= "GT_WU" times will be "
               "considered due to small memlimit.\n", maxfreq);
      }
    } else {
      gt_error_set(err,
                   "Option -memlimit too strict: need at least "GT_WU"MB",
                   (count >> 20) + 1);
      had_err = -1;
    }
    if (arg->verbose && !had_err) {
      gt_timer_show_formatted(vtimer,
                              "# Finished in "GT_WD".%06ld seconds. ",
                              stdout);
      if (maxfreq == GT_UWORD_MAX) {
        printf("Disable k-mer maximum frequency, ");
      } else {
        printf("Set k-mer maximum frequency to "GT_WU", ", maxfreq);
      }
      printf("expect "GT_WU" seed pairs.\n", count);
    }
  }

  /* create mlist of SeedPairs */
  if (!had_err) {
    if (arg->verbose) {
      printf("# Start building seed pairs on equal k-mers...\n");
      gt_timer_start(vtimer);
    }
    gt_diagbandseed_merge(&mlist,
                          alist,
                          alen,
                          blist,
                          blen,
                          &maxfreq,
                          maxgram,
                          arg->memlimit,
                          NULL, /* histogram not needed: save seed pairs */
                          endposdiff,
                          selfcomp);
  }

  gt_free(alist);
  if (!selfcomp || arg->mirror)
    gt_free(blist);

  /* sort mlist */
  if (!had_err) {
    mlen = mlist.nextfreeGtDiagbandseedSeedPair;
    rdxinfo = gt_radixsort_new_uint64keypair(mlen);
    gt_radixsort_inplace_Gtuint64keyPair((Gtuint64keyPair*)mlist.
                                         spaceGtDiagbandseedSeedPair,
                                         mlen);
    gt_radixsort_delete(rdxinfo);
    if (arg->verbose) {
      printf("# Collected and sorted "GT_WU" seed pairs ", mlen);
      gt_timer_show_formatted(vtimer, "in "GT_WD".%06ld seconds.\n", stdout);
    }
  }

  /* verify SeedPairs in the sequences */
  if ((arg->debug_seedpair || arg->verify) && mlen != 0) {
    GtDiagbandseedSeedPair *j = mlist.spaceGtDiagbandseedSeedPair;
    GtDiagbandseedSeedPair *last = j + mlen;
    char *buf1 = gt_malloc(3 * (1 + arg->seedlength) * sizeof *buf1);
    char *buf2 = buf1 + 1 + arg->seedlength;
    char *buf3 = buf2 + 1 + arg->seedlength;
    if (arg->verbose && arg->verify) {
      printf("# Start verifying seed pairs...\n");
      gt_timer_start(vtimer);
    }
    while (j < last && !had_err) {
      if (arg->debug_seedpair) {
        printf("SeedPair (%d,%d,%d,%d)\n",
               j->aseqnum, j->bseqnum, j->apos, j->bpos);
      }
      if (arg->verify) {
        char *idx = NULL;
        GtDiagbandseedPosition a = j->apos + gt_encseq_seqstartpos(aencseq,
                                                                   j->aseqnum);
        GtDiagbandseedPosition b = j->bpos + gt_encseq_seqstartpos(bencseq,
                                                                   j->bseqnum);
        gt_encseq_extract_decoded(aencseq, buf1, a + 1 - arg->seedlength, a);
        gt_encseq_extract_decoded(bencseq, buf2, b + 1 - arg->seedlength, b);
        buf1[arg->seedlength] = buf2[arg->seedlength] = '\0';
        buf3[arg->seedlength] = '\0';

        for (idx = buf3; idx < buf3 + arg->seedlength; idx++)
          gt_complement(idx, buf2[arg->seedlength + buf3 - idx - 1], NULL);
        if (strcmp(buf1, buf2) != 0 &&
            (!arg->mirror || strcmp(buf1, buf3) != 0))
        {
          gt_error_set(err, "Wrong seed(%d,%d,%d,%d): %s != %s / %s\n",
                       j->aseqnum, j->bseqnum, j->apos, j->bpos, buf1, buf2,
                       buf3);
          had_err = -1;
        }
      }
      j++;
    }
    gt_free(buf1);
    if (arg->verbose && !had_err && arg->verify) {
      printf("# Successfully verified each seed pair ");
      gt_timer_show_formatted(vtimer, "in "GT_WD".%06ld seconds.\n", stdout);
    }
  }

  /* process SeedPairs */
  if (had_err == 0 && arg->verbose &&
      (arg->extendgreedyinfo != NULL || arg->extendxdropinfo != NULL)) {
    printf("# Start seed pair extension...\n");
    gt_timer_start(vtimer);
  }
  if (had_err == 0 && mlen != 0) {
    had_err = gt_diagbandseed_process_seeds(aencseq,
                                            bencseq,
                                            &mlist,
                                            arg->extendgreedyinfo,
                                            arg->extendxdropinfo,
                                            arg->querymatchoutopt,
                                            arg->seedlength,
                                            arg->logdiagbandwidth,
                                            arg->mincoverage,
                                            amaxlen,
                                            bmaxlen,
                                            err);
  }
  GT_FREEARRAY(&mlist, GtDiagbandseedSeedPair);

  if (arg->verbose) {
    if (had_err == 0 &&
        (arg->extendgreedyinfo != NULL || arg->extendxdropinfo != NULL)) {
      gt_timer_show_formatted(vtimer, "# Completed extension in "GT_WD".%06ld "
                              "seconds.\n", stdout);
    }
    gt_timer_stop(vtimer);
    gt_timer_delete(vtimer);
  }

  return had_err;
}
