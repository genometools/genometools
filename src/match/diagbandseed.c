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

#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "core/codetype.h"
#include "core/complement.h"
#include "core/encseq.h"
#include "core/minmax.h"
#include "core/radix_sort.h"
#include "core/range_api.h"
#include "core/timer_api.h"
#include "core/warning_api.h"
#include "core/arraydef.h"
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

typedef struct {
  GtCodetype code;            /* only sort criterion */
  GtDiagbandseedSeqnum seqnum;
  GtDiagbandseedPosition endpos;
} GtDiagbandseedKmerPos;

typedef struct {
  GtDiagbandseedSeqnum bseqnum; /*  2nd important sort criterion */
  GtDiagbandseedSeqnum aseqnum; /* most important sort criterion */
  GtDiagbandseedPosition apos;
  GtDiagbandseedPosition bpos;  /*  3rd important sort criterion */
} GtDiagbandseedSeedPair;

GT_DECLAREARRAYSTRUCT(GtDiagbandseedSeedPair);

struct GtDiagbandseedProcKmerInfo {
  GtDiagbandseedKmerPos *list;
  GtUword numberofkmerscollected;
  GtDiagbandseedSeqnum seqnum;
  GtDiagbandseedPosition endpos;
  const GtEncseq *encseq;
  unsigned int seedlength;
  GtReadmode readmode;
  GtSpecialrangeiterator *sri;
  GtUword totallength;
  GtUword prev_separator;
  GtUword next_separator;
};

/* Returns the position of the next separator following separatorpos.
   If the end of encseq is reached, the position behind is returned. */
static GtUword gt_diagbandseed_update_separatorpos(GtUword separatorpos,
                                                   GtSpecialrangeiterator *sri,
                                                   const GtEncseq *encseq,
                                                   GtUword totallength,
                                                   GtReadmode readmode) {
  if (sri != NULL) {
    GtUword idx;
    GtRange range;

    gt_assert(encseq != NULL && separatorpos < totallength);
    range.start = separatorpos;
    range.end = separatorpos + 1;
    while (gt_specialrangeiterator_next(sri, &range)) {
      for (idx = range.start; idx < range.end; idx++) {
        if (gt_encseq_position_is_separator(encseq, idx, readmode)) {
          return idx;
        }
      }
    }
  }
  return totallength;
}

/* Add given code and its seqnum and position to a list. */
static void gt_diagbandseed_processkmercode(void *prockmerinfo,
                                            bool firstinrange,
                                            GtUword startpos,
                                            GtCodetype code)
{
  GtDiagbandseedProcKmerInfo *pkinfo;
  GtDiagbandseedKmerPos *kmerposptr;

  gt_assert(prockmerinfo != NULL);
  pkinfo = (GtDiagbandseedProcKmerInfo *) prockmerinfo;
  kmerposptr = pkinfo->list + pkinfo->numberofkmerscollected;

  /* check separator positions and determine next seqnum and endpos */
  if (firstinrange) {
    const GtUword endpos = startpos + pkinfo->seedlength - 1;
    while (endpos >= pkinfo->next_separator) {
      pkinfo->seqnum++;
      pkinfo->prev_separator = pkinfo->next_separator + 1;
      pkinfo->next_separator
        = gt_diagbandseed_update_separatorpos(pkinfo->next_separator,
                                              pkinfo->sri,
                                              pkinfo->encseq,
                                              pkinfo->totallength,
                                              pkinfo->readmode);
      gt_assert(pkinfo->next_separator > pkinfo->prev_separator);
    }
    gt_assert(endpos >= pkinfo->prev_separator);
    gt_assert(startpos < pkinfo->next_separator);
    if (pkinfo->readmode == GT_READMODE_FORWARD) {
      pkinfo->endpos = (GtDiagbandseedPosition) (endpos -
                                                 pkinfo->prev_separator);
    } else {
      pkinfo->endpos = (GtDiagbandseedPosition) (pkinfo->next_separator - 1 -
                                                 startpos);
    }
  }

  /* save k-mer code */
  kmerposptr->code
    = pkinfo->readmode == GT_READMODE_FORWARD
        ? code
        : gt_kmercode_reverse(code, pkinfo->seedlength);
  /* save endpos and seqnum */
  gt_assert(pkinfo->endpos != UINT_MAX);
  kmerposptr->endpos = pkinfo->endpos;
  pkinfo->endpos = (pkinfo->readmode == GT_READMODE_FORWARD
                    ? pkinfo->endpos + 1 : pkinfo->endpos - 1);
  kmerposptr->seqnum = pkinfo->seqnum;
  pkinfo->numberofkmerscollected++;
}

/* Uses GtKmercodeiterator for fetching the kmers. */
static void gt_diagbandseed_get_kmers_kciter(GtDiagbandseedProcKmerInfo *pkinfo)
{
  GtKmercodeiterator *kc_iter = NULL;
  const GtKmercode *kmercode = NULL;
  bool firstinrange = true;
  GtUword maxpos = 0, position;

  /* initialise GtKmercodeiterator */
  gt_assert(pkinfo != NULL);
  kc_iter = gt_kmercodeiterator_encseq_new(pkinfo->encseq,
                                           pkinfo->readmode,
                                           pkinfo->seedlength,
                                           0);
  if (pkinfo->seedlength <= pkinfo->totallength) {
    maxpos = pkinfo->totallength + 1 - pkinfo->seedlength;
  }

  /* iterate */
  for (position = 0; position < maxpos; position++)
  {
    kmercode = gt_kmercodeiterator_encseq_next(kc_iter);
    if (!kmercode->definedspecialposition) {
      gt_diagbandseed_processkmercode((void *) pkinfo,
                                      firstinrange,
                                      position,
                                      kmercode->code);
      firstinrange = false;
    } else {
      firstinrange = true;
    }
  }
  gt_kmercodeiterator_delete(kc_iter);
}

/* Returns a GtDiagbandseedKmerPos list of k-mers from a given encseq. */
static GtUword gt_diagbandseed_get_kmers(GtDiagbandseedKmerPos *list,
                                         const GtEncseq *encseq,
                                         unsigned int seedlength,
                                         GtReadmode readmode)
{
  GtDiagbandseedProcKmerInfo pkinfo;

  gt_assert(list != NULL && encseq != NULL);
  pkinfo.list = list;
  pkinfo.numberofkmerscollected = 0;
  pkinfo.seqnum = pkinfo.endpos = 0;
  pkinfo.encseq = encseq;
  pkinfo.seedlength = seedlength;
  pkinfo.readmode = readmode;
  pkinfo.totallength = gt_encseq_total_length(encseq);
  pkinfo.sri = NULL;
  if (gt_encseq_has_specialranges(encseq))
  {
    pkinfo.sri = gt_specialrangeiterator_new(encseq, true);
  }
  pkinfo.prev_separator = pkinfo.next_separator = 0;
  pkinfo.next_separator
    = gt_diagbandseed_update_separatorpos(pkinfo.next_separator,
                                          pkinfo.sri,
                                          encseq,
                                          pkinfo.totallength,
                                          readmode);
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
  } else
  {
    /* Use GtKmercodeiterator for encseq access */
    gt_diagbandseed_get_kmers_kciter(&pkinfo);
  }
  if (gt_encseq_has_specialranges(encseq))
  {
    gt_specialrangeiterator_delete(pkinfo.sri);
  }
  return pkinfo.numberofkmerscollected;
}

static void gt_diagbandseed_processhistogram(GtUword *histogram,
                                             GtUword *maxfreq,
                                             GtUword maxgram,
                                             GtUword memlimit,
                                             GtUword mem_used,
                                             bool alist_blist_id)
{
  /* calculate available memory, take 98% of memlimit */
  GtUword count = 0, frequency = 0;
  GtUword mem_avail = 0.98 * memlimit;
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
      gt_assert(frequency >= 2);
      frequency -= 2;
      gt_assert(count >= histogram[frequency]);
      count -= histogram[frequency];
    } else if (frequency == maxgram + 1) {
      frequency = GT_UWORD_MAX;
    }
    *maxfreq = MIN(*maxfreq, frequency);
  }

  /* determine minimum required memory for error message */
  if (*maxfreq <= 1 && alist_blist_id) {
    count = (histogram[0] + histogram[1]) * sizeof (GtDiagbandseedSeedPair);
    count = (count + mem_used) / 0.98;
  } else if (*maxfreq == 0) {
    count = histogram[0] * sizeof (GtDiagbandseedSeedPair);
    count = (count + mem_used) / 0.98;
  }
  histogram[maxgram] = count;
}

/* Returns a GtDiagbandseedSeedPair list of equal kmers from lists a and b. */
static void gt_diagbandseed_merge(GtArrayGtDiagbandseedSeedPair *mlist,
                                  const GtDiagbandseedKmerPos *alist,
                                  GtUword alen,
                                  const GtDiagbandseedKmerPos *blist,
                                  GtUword blen,
                                  GtUword *maxfreq,
                                  GtUword maxgram,
                                  GtUword memlimit,
                                  GtUword *histogram,
                                  unsigned int endposdiff,
                                  bool selfcomp,
                                  GtUword len_used)
{
  const GtDiagbandseedKmerPos *aptr = alist, *bptr = blist, *aend, *bend;
  const GtUword array_incr = 256;
  GtUword frequency = 0;

  gt_assert(alist != NULL && blist != NULL && maxfreq != NULL);
  gt_assert((histogram == NULL && mlist != NULL) ||
            (histogram != NULL && mlist == NULL));
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
    gt_diagbandseed_processhistogram(histogram,
                                     maxfreq,
                                     maxgram,
                                     memlimit,
                                     len_used * sizeof *alist,
                                     alist == blist);
  }
}

/* Verify seed pairs in the original sequences */
int gt_diagbandseed_verify(const GtEncseq *aencseq,
                           const GtEncseq *bencseq,
                           const GtArrayGtDiagbandseedSeedPair *mlist,
                           unsigned int seedlength,
                           bool reverse,
                           GtError *err) {
  GtDiagbandseedSeedPair *curr_sp, *max_sp;
  char *buf1 = gt_malloc(3 * (seedlength + 1) * sizeof *buf1);
  char *buf2 = buf1 + 1 + seedlength;
  char *buf3 = buf2 + 1 + seedlength;
  buf1[seedlength] = buf2[seedlength] = buf3[seedlength] = '\0';

  gt_assert(mlist != NULL && aencseq != NULL && bencseq != NULL);
  curr_sp = mlist->spaceGtDiagbandseedSeedPair;
  max_sp = curr_sp + mlist->nextfreeGtDiagbandseedSeedPair;

  while (curr_sp < max_sp) {
    GtDiagbandseedPosition apos, bpos;

    /* extract decoded k-mers at seed pair positions */
    apos = curr_sp->apos + gt_encseq_seqstartpos(aencseq, curr_sp->aseqnum);
    gt_encseq_extract_decoded(aencseq, buf1, apos + 1 - seedlength, apos);

    if (!reverse) {
      bpos = curr_sp->bpos + gt_encseq_seqstartpos(bencseq, curr_sp->bseqnum);
      gt_encseq_extract_decoded(bencseq, buf2, bpos + 1 - seedlength, bpos);
      if (strcmp(buf1, buf2) != 0) {
        gt_error_set(err, "Wrong SeedPair (%d,%d,%d,%d): %s != %s\n",
                     curr_sp->aseqnum, curr_sp->bseqnum, curr_sp->apos,
                     curr_sp->bpos, buf1, buf2);
        gt_free(buf1);
        return -1;
      }
    } else {
      /* get reverse k-mer */
      char *idx;
      bpos = gt_encseq_seqstartpos(bencseq, curr_sp->bseqnum) +
             gt_encseq_seqlength(bencseq, curr_sp->bseqnum) - curr_sp->bpos - 1;
      gt_encseq_extract_decoded(bencseq, buf2, bpos, bpos + seedlength - 1);

      for (idx = buf3; idx < buf3 + seedlength; idx++) {
        gt_complement(idx, buf2[seedlength + buf3 - idx - 1], NULL);
      }
      if (strcmp(buf1, buf3) != 0) {
        gt_error_set(err, "Wrong SeedPair (%d,%d,%d,%d): %s != %s\n",
                     curr_sp->aseqnum, curr_sp->bseqnum, curr_sp->apos,
                     curr_sp->bpos, buf1, buf3);
        gt_free(buf1);
        return -1;
      }
    }
    curr_sp++;
  }
  gt_free(buf1);
  return 0;
}

/* start seed extension for seed pairs in mlist */
static GtUword gt_diagbandseed_process_seeds(const GtEncseq *aencseq,
                                  const GtEncseq *bencseq,
                                  const GtArrayGtDiagbandseedSeedPair *mlist,
                                  GtGreedyextendmatchinfo *extendgreedyinfo,
                                  GtXdropmatchinfo *extendxdropinfo,
                                  GtQuerymatchoutoptions *querymatchoutopt,
                                  bool seed_display,
                                  unsigned int seedlength,
                                  GtUword errorpercentage,
                                  GtUword userdefinedleastlength,
                                  GtUword logdiagbandwidth,
                                  GtUword mincoverage,
                                  GtUword amaxlen,
                                  GtUword bmaxlen,
                                  bool reverse)
{
  GtDiagbandseedScore *score = NULL;
  GtDiagbandseedPosition *lastp = NULL;
  GtExtendSelfmatchRelativeFunc extend_selfmatch_relative_function = NULL;
  GtExtendQuerymatchRelativeFunc extend_querymatch_relative_function = NULL;
  GtProcessinfo_and_querymatchspaceptr info_querymatch;
  const GtDiagbandseedSeedPair *lm = NULL, *maxsegm, *nextsegm, *idx;
  const GtUword ndiags = (amaxlen >> logdiagbandwidth) +
                         (bmaxlen >> logdiagbandwidth) + 2;
  const GtUword minsegmentlen = (mincoverage - 1) / seedlength + 1;
  GtUword mlen = 0, diag = 0;
  bool firstinrange = true;
  GtUword count_extensions = 0;
  const GtReadmode query_readmode
    = (reverse && (extendgreedyinfo != NULL ||
                   extendxdropinfo != NULL)) ? GT_READMODE_REVCOMPL
                                             : GT_READMODE_FORWARD;

  gt_assert(mlist != NULL);
  mlen = mlist->nextfreeGtDiagbandseedSeedPair; /* mlist length  */
  lm = mlist->spaceGtDiagbandseedSeedPair;      /* mlist pointer */

  /* select extension method */
  if (extendgreedyinfo != NULL) {
    info_querymatch.processinfo = (void *) extendgreedyinfo;
    extend_selfmatch_relative_function = gt_greedy_extend_selfmatch_relative;
    extend_querymatch_relative_function = gt_greedy_extend_querymatch_relative;
  } else if (extendxdropinfo != NULL) {
    info_querymatch.processinfo = (void *) extendxdropinfo;
    extend_selfmatch_relative_function = gt_xdrop_extend_selfmatch_relative;
    extend_querymatch_relative_function = gt_xdrop_extend_querymatch_relative;
  } else { /* no seed extension */
    return 0;
  }

  if (mlen < minsegmentlen)
    return 0;
  gt_assert(aencseq != NULL && bencseq != NULL);

  info_querymatch.querymatchspaceptr = gt_querymatch_new();
  if (seed_display)
  {
    gt_querymatch_seed_display_set(info_querymatch.querymatchspaceptr);
  }
  if (querymatchoutopt != NULL)
  {
    gt_querymatch_outoptions_set(info_querymatch.querymatchspaceptr,
                                 querymatchoutopt);
  }
  gt_querymatch_query_readmode_set(info_querymatch.querymatchspaceptr,
                                   query_readmode);

  /* score[0] and score[ndiags+1] remain zero for boundary */
  score = gt_calloc(ndiags + 2, sizeof *score);
  lastp = gt_calloc(ndiags, sizeof *lastp);
  maxsegm = lm + mlen - minsegmentlen;
  nextsegm = lm;

  /* iterate through segments of equal k-mers */
  while (nextsegm <= maxsegm) {
    const GtDiagbandseedSeedPair *currsegm = nextsegm;

    /* if insuffienct number of kmers in segment: skip whole segment */
    if (currsegm->aseqnum != (currsegm + minsegmentlen - 1)->aseqnum ||
        currsegm->bseqnum != (currsegm + minsegmentlen - 1)->bseqnum) {
      do {
        nextsegm++;
      } while (nextsegm < lm + mlen &&
               nextsegm->aseqnum == currsegm->aseqnum &&
               nextsegm->bseqnum == currsegm->bseqnum);
      continue;
    }

    /* calculate diagonal band scores */
    do {
      gt_assert(nextsegm->bpos <= bmaxlen && nextsegm->apos <= amaxlen);
      diag = (amaxlen + (GtUword)nextsegm->bpos - (GtUword)nextsegm->apos)
             >> logdiagbandwidth;
      if (nextsegm->bpos >= seedlength + lastp[diag]) {
        /* no overlap: add seedlength */
        score[diag + 1] += seedlength;
      } else {
        /* overlap: add difference below overlap */
        gt_assert(lastp[diag] <= nextsegm->bpos); /* if fail: sort by bpos */
        score[diag + 1] += nextsegm->bpos - lastp[diag];
      }
      lastp[diag] = nextsegm->bpos;
      nextsegm++;
    } while (nextsegm < lm + mlen && nextsegm->aseqnum == currsegm->aseqnum &&
             nextsegm->bseqnum == currsegm->bseqnum);

    /* test for mincoverage and overlap to previous extension */
    firstinrange = true;
    for (idx = currsegm; idx < nextsegm; idx++) {
      gt_assert(idx->apos <= amaxlen);
      diag = (amaxlen + (GtUword)idx->bpos - (GtUword)idx->apos)
             >> logdiagbandwidth;
      if ((GtUword)MAX(score[diag + 2], score[diag]) + (GtUword)score[diag + 1]
          >= mincoverage)
      {
        /* relative seed start position in B */
        const GtUword bstart = (GtUword) (idx->bpos + 1 - seedlength);

        if (firstinrange ||
            !gt_querymatch_overlap(info_querymatch.querymatchspaceptr,bstart))
        {
          /* relative seed start position in A */
          const GtUword astart = (GtUword) (idx->apos + 1 - seedlength);
          /* extend seed */
          const GtQuerymatch *querymatch = NULL;
          if (aencseq == bencseq) {
            count_extensions++;
            querymatch = extend_selfmatch_relative_function(&info_querymatch,
                                                            aencseq,
                                                            idx->aseqnum,
                                                            astart,
                                                            idx->bseqnum,
                                                            bstart,
                                                            seedlength,
                                                            query_readmode);
          } else {
            count_extensions++;
            querymatch = extend_querymatch_relative_function(&info_querymatch,
                                                             aencseq,
                                                             idx->aseqnum,
                                                             astart,
                                                             bencseq,
                                                             idx->bseqnum,
                                                             bstart,
                                                             seedlength,
                                                             query_readmode);
          }
          if (querymatch != NULL) {
            firstinrange = false;
            /* show extension results */
            if (gt_querymatch_check_final(querymatch,errorpercentage,
                                          userdefinedleastlength))
            {
              gt_querymatch_prettyprint(querymatch);
            }
          }
        }
      }
    }

    /* reset diagonal band scores */
    for (idx = currsegm; idx < nextsegm; idx++) {
      diag = (amaxlen + (GtUword)idx->bpos - (GtUword)idx->apos)
             >> logdiagbandwidth;
      score[diag + 1] = 0;
      lastp[diag] = 0;
    }
  }
  gt_querymatch_delete(info_querymatch.querymatchspaceptr);
  gt_free(score);
  gt_free(lastp);
  return count_extensions;
}

static GtUword gt_seed_extend_numofkmers(const GtEncseq *encseq,
                                         GtUword seedlength)
{
  GtUword kmers;
  const GtUword totallength = gt_encseq_total_length(encseq),
                min_seq_length = gt_encseq_min_seq_length(encseq),
                num_of_sequences = gt_encseq_num_of_sequences(encseq),
                num_specialchar = gt_encseq_specialcharacters(encseq);
  gt_assert(seedlength > 0);
  kmers = totallength -
          (num_of_sequences * MIN(seedlength, min_seq_length + 1) - 1);
  kmers = MIN(kmers, totallength - (seedlength - 1 + num_specialchar));
  return kmers;
}

static int gt_diagbandseed_prepare_mlist2(GtArrayGtDiagbandseedSeedPair *mlist,
                                          GtDiagbandseedKmerPos *alist,
                                          GtUword alen,
                                          GtUword blen,
                                          const GtEncseq *aencseq,
                                          const GtEncseq *bencseq,
                                          const GtDiagbandseed *arg,
                                          GtUword maxfreq,
                                          bool selfcomp,
                                          GtError *err)
{
  GtDiagbandseedKmerPos *clist;
  GtRadixsortinfo *rdxinfo;
  GtTimer *vtimer = NULL;
  GtUword mlen, clen;
  int had_err = 0;

  clen = selfcomp ? alen : blen;
  if (arg->verbose) {
    vtimer = gt_timer_new();
    printf("# Start fetching " GT_WU " reverse complement %u-mers...\n",
           clen, arg->seedlength);
    gt_timer_start(vtimer);
  }

  clist = gt_malloc(clen * sizeof *clist);
  clen = gt_diagbandseed_get_kmers(clist,
                                   selfcomp ? aencseq : bencseq,
                                   arg->seedlength,
                                   GT_READMODE_COMPL);
  if (arg->debug_kmer) {
    GtDiagbandseedKmerPos *idx;
    for (idx = clist; idx < clist + clen; idx++) {
      printf("# Kmer (" GT_LX ",%d,%d)\n",
             idx->code, idx->endpos, idx->seqnum);
    }
  }

  if (arg->verbose) {
    printf("# ...found " GT_WU " %u-mers ", clen, arg->seedlength);
    gt_timer_show_formatted(vtimer, "in " GT_WD ".%06ld seconds.\n", stdout);
    gt_timer_start(vtimer);
  }

  rdxinfo = gt_radixsort_new_ulongpair(clen);
  gt_radixsort_inplace_GtUwordPair((GtUwordPair *)clist, clen);
  gt_radixsort_delete(rdxinfo);

  if (arg->verbose) {
    printf("# ...sorted " GT_WU " %u-mers ", clen, arg->seedlength);
    gt_timer_show_formatted(vtimer, "in " GT_WD ".%06ld seconds.\n", stdout);
  }

  if (arg->verbose) {
    printf("# Start building seed pairs using rev.compl. %u-mers...\n",
           arg->seedlength);
    gt_timer_start(vtimer);
  }

  GT_INITARRAY(mlist, GtDiagbandseedSeedPair);
  gt_diagbandseed_merge(mlist,
                        alist,
                        alen,
                        clist,
                        clen,
                        &maxfreq,
                        1,    /* maxgram not needed */
                        arg->memlimit,
                        NULL, /* histogram not needed: save seed pairs */
                        0,
                        selfcomp,
                        0);   /* len_used not needed */
  mlen = mlist->nextfreeGtDiagbandseedSeedPair;
  gt_free(alist);
  alist = NULL;
  gt_free(clist);
  clist = NULL;

  if (arg->verbose) {
    printf("# ...collected " GT_WU " rev.compl. seed pairs ", mlen);
    gt_timer_show_formatted(vtimer, "in " GT_WD ".%06ld seconds.\n", stdout);
  }

  if (mlen > 0) {
    if (arg->verbose) {
      gt_timer_start(vtimer);
    }
    rdxinfo = gt_radixsort_new_uint64keypair(mlen);
    gt_radixsort_inplace_Gtuint64keyPair((Gtuint64keyPair*) mlist->
                                         spaceGtDiagbandseedSeedPair,
                                         mlen);
    gt_radixsort_delete(rdxinfo);
    if (arg->verbose) {
      printf("# ...sorted " GT_WU " rev.compl. seed pairs in ", mlen);
      gt_timer_show_formatted(vtimer, GT_WD ".%06ld seconds.\n", stdout);
    }

    if (arg->debug_seedpair) {
      GtDiagbandseedSeedPair *curr_sp = mlist->spaceGtDiagbandseedSeedPair;
      while (curr_sp < mlist->spaceGtDiagbandseedSeedPair + mlen) {
        printf("# SeedPair (%d,%d,%d,%d)\n", curr_sp->aseqnum,
               curr_sp->bseqnum, curr_sp->apos, curr_sp->bpos);
        curr_sp++;
      }
    }
  }

  if (arg->verify) {
    if (arg->verbose) {
      printf("# Start verifying rev.compl. seed pairs...\n");
      gt_timer_start(vtimer);
    }

    had_err = gt_diagbandseed_verify(aencseq,
                                     bencseq,
                                     mlist,
                                     arg->seedlength,
                                     true,
                                     err);
    if (had_err) {
      GT_FREEARRAY(mlist, GtDiagbandseedSeedPair);
    }

    if (!had_err && arg->verbose) {
      printf("# ...successfully verified each seed pair ");
      gt_timer_show_formatted(vtimer, "in " GT_WD ".%06ld seconds.\n", stdout);
    }
  }
  if (arg->verbose) {
    gt_timer_delete(vtimer);
  }
  return had_err;
}

int gt_diagbandseed_run(const GtEncseq *aencseq,
                        const GtEncseq *bencseq,
                        const GtDiagbandseed *arg,
                        GtError *err)
{
  GtDiagbandseedKmerPos *alist = NULL, *blist = NULL;
  GtArrayGtDiagbandseedSeedPair mlist, mlist_rev;
  GtTimer *vtimer = NULL;
  GtRadixsortinfo *rdxinfo;
  GtUword alen, blen = 0, mlen = 0, maxfreq, amaxlen, bmaxlen, ankmers;
  unsigned int endposdiff;
  int had_err = 0;
  bool alist_blist_id = false;
  bool both_strands = false;
  const GtUword maxgram = 10000; /* Cap on k-mer count histogram */
  const bool selfcomp = (bencseq == aencseq) ? true : false;

  gt_assert(arg != NULL && aencseq != NULL && bencseq != NULL);
  maxfreq = arg->maxfreq;
  endposdiff = !arg->overlappingseeds ? arg->seedlength : 1;
  amaxlen = gt_encseq_max_seq_length(aencseq);
  bmaxlen = gt_encseq_max_seq_length(bencseq);
  alist_blist_id = (selfcomp && !arg->nofwd) ? true : false;
  both_strands = (arg->norev || arg->nofwd) ? false : true;
  if (MIN(amaxlen, bmaxlen) < arg->seedlength) {
    gt_error_set(err,"argument to option \"-seedlength\" must be an integer <= "
                 GT_WU " (length of longest sequence).", MIN(amaxlen, bmaxlen));
    return -1;
  }

  /* estimate number of expected kmers */
  ankmers = gt_seed_extend_numofkmers(aencseq, arg->seedlength);

  /* prepare list of kmers from aencseq */
  if (arg->verbose) {
    vtimer = gt_timer_new();
    printf("# Start fetching (at most " GT_WU ") %u-mers for list A...\n",
           ankmers,arg->seedlength);
    gt_timer_start(vtimer);
  }

  alist = gt_malloc(ankmers * sizeof *alist);
  alen = gt_diagbandseed_get_kmers(alist,
                                   aencseq,
                                   arg->seedlength,
                                   GT_READMODE_FORWARD);
  if (arg->debug_kmer) {
    GtDiagbandseedKmerPos *idx;
    for (idx = alist; idx < alist + alen; idx++) {
      printf("# Kmer (" GT_LX ",%d,%d)\n", idx->code, idx->endpos, idx->seqnum);
    }
  }

  if (arg->verbose) {
    printf("# ...found " GT_WU " %u-mers ", alen, arg->seedlength);
    gt_timer_show_formatted(vtimer, "in " GT_WD ".%06ld seconds.\n", stdout);
    gt_timer_start(vtimer);
  }

  /* sort alist */
  rdxinfo = gt_radixsort_new_ulongpair(alen);
  gt_radixsort_inplace_GtUwordPair((GtUwordPair *)alist, alen);
  gt_radixsort_delete(rdxinfo);

  if (arg->verbose) {
    printf("# ...sorted " GT_WU " %u-mers ", alen, arg->seedlength);
    gt_timer_show_formatted(vtimer, "in " GT_WD ".%06ld seconds.\n", stdout);
  }

  /* if necessary: prepare list of kmers from bencseq */
  if (alist_blist_id) {
    /* compare reads of encseq A with themselves */
    blist = alist;
    blen = alen;
  } else {
    GtRadixsortinfo *rdxinfo;
    const GtReadmode readmode = arg->nofwd ? GT_READMODE_COMPL
                                           : GT_READMODE_FORWARD;
    GtUword bnkmers = 0;

    /* estimate number of expected kmers */
    if (selfcomp) {                /* same kmer amount as for alist */
      bnkmers = alen;
    } else {                 /* no selfcomp: calculate from bencseq */
      bnkmers = gt_seed_extend_numofkmers(bencseq, arg->seedlength);
    }

    if (arg->verbose) {
      printf("# Start fetching (at most " GT_WU ") %u-mers for list B...\n",
             bnkmers,arg->seedlength);
      gt_timer_start(vtimer);
    }

    /* fill blist */
    blist = gt_malloc(bnkmers * sizeof *blist);
    blen = gt_diagbandseed_get_kmers(blist,
                                     selfcomp ? aencseq : bencseq,
                                     arg->seedlength,
                                     readmode);
    if (arg->debug_kmer) {
      GtDiagbandseedKmerPos *idx;
      for (idx = blist; idx < blist + blen; idx++) {
        printf("# Kmer (" GT_LX ",%d,%d)\n", idx->code, idx->endpos,
               idx->seqnum);
      }
    }

    if (arg->verbose) {
      printf("# ...found " GT_WU " %u-mers ", blen, arg->seedlength);
      gt_timer_show_formatted(vtimer, "in " GT_WD ".%06ld seconds.\n", stdout);
      gt_timer_start(vtimer);
    }

    /* sort blist */
    rdxinfo = gt_radixsort_new_ulongpair(blen);
    gt_radixsort_inplace_GtUwordPair((GtUwordPair *)blist, blen);
    gt_radixsort_delete(rdxinfo);

    if (arg->verbose) {
      printf("# ...sorted " GT_WU " %u-mers ", blen, arg->seedlength);
      gt_timer_show_formatted(vtimer, "in " GT_WD ".%06ld seconds.\n", stdout);
    }
  }

  /* calculate maxfreq from memlimit */
  if (arg->memlimit < GT_UWORD_MAX) {
    GtUword *histogram = NULL;
    GtUword len_used = alen;
    if (!selfcomp || !arg->norev) {
      len_used += blen;
    }
    if (!selfcomp && both_strands && arg->extend_last) {
      len_used += blen;
    }

    if (arg->verbose) {
      printf("# Start calculating k-mer frequency histogram...\n");
      gt_timer_start(vtimer);
    }

    /* build histogram; histogram[maxgram] := estimation for mlen */
    histogram = gt_calloc(maxgram + 1, sizeof *histogram);
    gt_diagbandseed_merge(NULL, /* mlist not needed: just count */
                          alist,
                          alen,
                          blist,
                          blen,
                          &maxfreq,
                          maxgram,
                          arg->memlimit,
                          histogram,
                          alist_blist_id ? endposdiff : 0,
                          selfcomp,
                          len_used);
    mlen = histogram[maxgram];
    gt_free(histogram);

    /* check maxfreq value */
    if (maxfreq == 0 || (maxfreq == 1 && alist_blist_id)) {
      gt_error_set(err,
                   "option -memlimit too strict: need at least " GT_WU "MB",
                   (mlen >> 20) + 1);
      mlen = 0;
      had_err = -1;
    } else if (maxfreq < 10) {
      gt_warning("Only %u-mers occurring <= " GT_WU " times will be "
                 "considered, "
                 "due to small memlimit. Expect " GT_WU " seed pairs.",
                 arg->seedlength,maxfreq, mlen);
    } else if (arg->verbose) {
      if (maxfreq == GT_UWORD_MAX) {
        printf("# Disable k-mer maximum frequency, ");
      } else {
        printf("# Set k-mer maximum frequency to " GT_WU ", ", maxfreq);
      }
      printf("expect " GT_WU " seed pairs.\n", mlen);
    }

    if (!had_err && arg->verbose) {
      gt_timer_show_formatted(vtimer,
                              "# ...finished in " GT_WD ".%06ld seconds.\n",
                              stdout);
    }
  } /* end of maxfreq calculation */

  /* create mlist of SeedPairs by merging alist and blist */
  if (!had_err) {
    GT_INITARRAY(&mlist, GtDiagbandseedSeedPair);
    if (mlen > 0) {
      /* allocate mlist space according to seed pair count */
      GT_CHECKARRAYSPACEMULTI(&mlist, GtDiagbandseedSeedPair, mlen);
    }

    if (arg->verbose) {
      printf("# Start building seed pairs on equal %u-mers...\n",
              arg->seedlength);
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
                          alist_blist_id ? endposdiff : 0,
                          selfcomp,
                          0);   /* len_used not needed */
    mlen = mlist.nextfreeGtDiagbandseedSeedPair;
    if (!both_strands) {
      gt_free(alist);
      alist = NULL;
    }
    if (!alist_blist_id) {
      gt_free(blist);
      blist = NULL;
    }

    if (arg->verbose) {
      printf("# ...collected " GT_WU " seed pairs ", mlen);
      gt_timer_show_formatted(vtimer, "in " GT_WD ".%06ld seconds.\n", stdout);
    }

    if (mlen > 0) {
      GtRadixsortinfo *rdxinfo;
      if (arg->verbose) {
        gt_timer_start(vtimer);
      }
      rdxinfo = gt_radixsort_new_uint64keypair(mlen);
      gt_radixsort_inplace_Gtuint64keyPair((Gtuint64keyPair*) mlist.
                                           spaceGtDiagbandseedSeedPair,
                                           mlen);
      gt_radixsort_delete(rdxinfo);
      if (arg->verbose) {
        printf("# ...sorted " GT_WU " seed pairs in ", mlen);
        gt_timer_show_formatted(vtimer, GT_WD ".%06ld seconds.\n", stdout);
      }

      if (arg->debug_seedpair) {
        GtDiagbandseedSeedPair *curr_sp = mlist.spaceGtDiagbandseedSeedPair;
        while (curr_sp < mlist.spaceGtDiagbandseedSeedPair + mlen) {
          printf("# SeedPair (%d,%d,%d,%d)\n", curr_sp->aseqnum,
                 curr_sp->bseqnum, curr_sp->apos, curr_sp->bpos);
          curr_sp++;
        }
      }
    }

    /* verify SeedPairs in the sequences */
    if (!had_err && arg->verify) {
      if (arg->verbose) {
        printf("# Start verifying seed pairs...\n");
        gt_timer_start(vtimer);
      }

      had_err = gt_diagbandseed_verify(aencseq,
                                       bencseq,
                                       &mlist,
                                       arg->seedlength,
                                       arg->nofwd,
                                       err);
      if (had_err) {
        GT_FREEARRAY(&mlist, GtDiagbandseedSeedPair);
      }

      if (!had_err && arg->verbose) {
        printf("# ...successfully verified each seed pair in ");
        gt_timer_show_formatted(vtimer, GT_WD ".%06ld seconds.\n", stdout);
      }
    }
  } else { /* had_err */
    gt_free(alist);
    alist = NULL;
    if (!alist_blist_id) {
      gt_free(blist);
      blist = NULL;
    }
  }

  if (!had_err && both_strands && arg->extend_last) {
    had_err = gt_diagbandseed_prepare_mlist2(&mlist_rev,
                                             alist,
                                             alen,
                                             blen,
                                             aencseq,
                                             bencseq,
                                             arg,
                                             maxfreq,
                                             selfcomp,
                                             err);
  }

  /* process SeedPairs */
  if (!had_err && mlen > 0) {
    GtUword count_seedextensions;
    if (arg->verbose &&
        (arg->extendgreedyinfo != NULL || arg->extendxdropinfo != NULL)) {
      printf("# Start seed pair extension...\n");
      gt_timer_start(vtimer);
    }

    count_seedextensions
      = gt_diagbandseed_process_seeds(aencseq,
                                      bencseq,
                                      &mlist,
                                      arg->extendgreedyinfo,
                                      arg->extendxdropinfo,
                                      arg->querymatchoutopt,
                                      arg->seed_display,
                                      arg->seedlength,
                                      arg->errorpercentage,
                                      arg->userdefinedleastlength,
                                      arg->logdiagbandwidth,
                                      arg->mincoverage,
                                      amaxlen,
                                      bmaxlen,
                                      arg->nofwd);
    GT_FREEARRAY(&mlist, GtDiagbandseedSeedPair);
    if (!had_err && arg->verbose &&
        (arg->extendgreedyinfo != NULL || arg->extendxdropinfo != NULL)) {
      printf("# ...finished " GT_WU " seed pair extension%s ",
             count_seedextensions,
             count_seedextensions > 1 ? "s" : "");
      gt_timer_show_formatted(vtimer, "in " GT_WD ".%06ld seconds.\n", stdout);
    }
  }

  if (!had_err && both_strands && !arg->extend_last) {
    had_err = gt_diagbandseed_prepare_mlist2(&mlist_rev,
                                             alist,
                                             alen,
                                             blen,
                                             aencseq,
                                             bencseq,
                                             arg,
                                             maxfreq,
                                             selfcomp,
                                             err);
  }

  if (!had_err && both_strands && mlist_rev.nextfreeGtDiagbandseedSeedPair > 0)
  {
    GtUword count_seedextensions;
    if (arg->verbose &&
        (arg->extendgreedyinfo != NULL || arg->extendxdropinfo != NULL)) {
      gt_timer_start(vtimer);
    }
    count_seedextensions
      = gt_diagbandseed_process_seeds(aencseq,
                                      bencseq,
                                      &mlist_rev,
                                      arg->extendgreedyinfo,
                                      arg->extendxdropinfo,
                                      arg->querymatchoutopt,
                                      arg->seed_display,
                                      arg->seedlength,
                                      arg->errorpercentage,
                                      arg->userdefinedleastlength,
                                      arg->logdiagbandwidth,
                                      arg->mincoverage,
                                      amaxlen,
                                      bmaxlen,
                                      true);
    GT_FREEARRAY(&mlist_rev, GtDiagbandseedSeedPair);
    if (!had_err && arg->verbose &&
        (arg->extendgreedyinfo != NULL || arg->extendxdropinfo != NULL)) {
      printf("# ...finished " GT_WU " rev.compl. seed pair extension%s in ",
             count_seedextensions,
             count_seedextensions > 1 ? "s" : "");
      gt_timer_show_formatted(vtimer, GT_WD ".%06ld seconds.\n", stdout);
    }
  }

  if (arg->verbose) {
    gt_timer_delete(vtimer);
  }
  return had_err;
}
