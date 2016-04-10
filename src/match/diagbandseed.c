/*
  Copyright (c) 2015-2016 Joerg Winkler <joerg.winkler@studium.uni-hamburg.de>
  Copyright (c) 2015-2016 Center for Bioinformatics, University of Hamburg

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
#include "core/arraydef.h"
#include "core/codetype.h"
#include "core/complement.h"
#include "core/encseq.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/radix_sort.h"
#include "core/timer_api.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "match/diagbandseed.h"
#include "match/ft-front-prune.h"
#include "match/kmercodes.h"
#include "match/querymatch.h"
#include "match/querymatch-align.h"
#include "match/seed-extend.h"
#include "match/sfx-mappedstr.h"
#include "match/sfx-suffixer.h"

#ifdef GT_THREADS_ENABLED
#include "core/multithread_api.h"
#endif

#define GT_DIAGBANDSEED_SEQNUM_UNDEF UINT_MAX
/* #define GT_DIAGBANDSEED_SEEDHISTOGRAM 100 */

typedef uint32_t GtDiagbandseedPosition;
typedef uint32_t GtDiagbandseedSeqnum;
typedef uint32_t GtDiagbandseedScore;
typedef struct GtDiagbandseedProcKmerInfo GtDiagbandseedProcKmerInfo;

typedef struct {
  GtCodetype code;              /* only sort criterion */
  GtDiagbandseedSeqnum seqnum;
  GtDiagbandseedPosition endpos;
} GtDiagbandseedKmerPos;

GT_DECLAREARRAYSTRUCT(GtDiagbandseedKmerPos);

typedef struct {
  GtDiagbandseedSeqnum bseqnum; /*  2nd important sort criterion */
  GtDiagbandseedSeqnum aseqnum; /* most important sort criterion */
  GtDiagbandseedPosition apos;
  GtDiagbandseedPosition bpos;  /*  3rd important sort criterion */
} GtDiagbandseedSeedPair;

GT_DECLAREARRAYSTRUCT(GtDiagbandseedSeedPair);

struct GtDiagbandseedInfo {
  GtEncseq *aencseq;
  GtEncseq *bencseq;
  GtUword maxfreq;
  GtUword memlimit;
  unsigned int seedlength;
  bool norev;
  bool nofwd;
  bool overlappingseeds;
  bool verify;
  bool verbose;
  bool debug_kmer;
  bool debug_seedpair;
  bool extend_last;
  GtDiagbandseedExtendParams *extp;
};

struct GtDiagbandseedExtendParams {
  GtUword errorpercentage;
  GtUword userdefinedleastlength;
  GtUword logdiagbandwidth;
  GtUword mincoverage;
  unsigned int display_flag;
  bool use_apos;
  GtXdropscore xdropbelowscore;
  bool extendgreedy;
  bool extendxdrop;
  GtUword maxalignedlendifference;
  GtUword history_size;
  GtUword perc_mat_history;
  GtExtendCharAccess extend_char_access;
  GtUword sensitivity;
  double matchscore_bias;
  bool weakends;
  bool benchmark;
  GtUword alignmentwidth;
  bool always_polished_ends;
};

struct GtDiagbandseedProcKmerInfo {
  GtArrayGtDiagbandseedKmerPos *list;
  GtDiagbandseedSeqnum seqnum;
  GtDiagbandseedPosition endpos;
  const GtEncseq *encseq;
  unsigned int seedlength;
  GtReadmode readmode;
  GtSpecialrangeiterator *sri;
  GtUword totallength;
  GtUword prev_separator;
  GtUword next_separator;
  GtRange *specialrange;
};

/* * * * * CONSTRUCTORS AND DESTRUCTORS * * * * */

GtDiagbandseedInfo *gt_diagbandseed_info_new(GtEncseq *aencseq,
                                             GtEncseq *bencseq,
                                             GtUword maxfreq,
                                             GtUword memlimit,
                                             unsigned int seedlength,
                                             bool norev,
                                             bool nofwd,
                                             bool overlappingseeds,
                                             bool verify,
                                             bool verbose,
                                             bool debug_kmer,
                                             bool debug_seedpair,
                                             bool extend_last,
                                             GtDiagbandseedExtendParams *extp)
{
  GtDiagbandseedInfo *info = gt_malloc(sizeof *info);
  info->aencseq = aencseq;
  info->bencseq = bencseq;
  info->maxfreq = maxfreq;
  info->memlimit = memlimit;
  info->seedlength = seedlength;
  info->norev = norev;
  info->nofwd = nofwd;
  info->overlappingseeds = overlappingseeds;
  info->verify = verify;
  info->verbose = verbose;
  info->debug_kmer = debug_kmer;
  info->debug_seedpair = debug_seedpair;
  info->extend_last = extend_last;
  info->extp = extp;
  return info;
}

void gt_diagbandseed_info_delete(GtDiagbandseedInfo *info)
{
  if (info != NULL) {
    gt_encseq_delete(info->aencseq);
    gt_encseq_delete(info->bencseq);
    gt_diagbandseed_extend_params_delete(info->extp);
    gt_free(info);
  }
}

GtDiagbandseedExtendParams *gt_diagbandseed_extend_params_new(
                                GtUword errorpercentage,
                                GtUword userdefinedleastlength,
                                GtUword logdiagbandwidth,
                                GtUword mincoverage,
                                unsigned int display_flag,
                                bool use_apos,
                                GtXdropscore xdropbelowscore,
                                bool extendgreedy,
                                bool extendxdrop,
                                GtUword maxalignedlendifference,
                                GtUword history_size,
                                GtUword perc_mat_history,
                                GtExtendCharAccess extend_char_access,
                                GtUword sensitivity,
                                double matchscore_bias,
                                bool weakends,
                                bool benchmark,
                                GtUword alignmentwidth,
                                bool always_polished_ends)
{
  GtDiagbandseedExtendParams *extp = gt_malloc(sizeof *extp);
  extp->errorpercentage = errorpercentage;
  extp->userdefinedleastlength = userdefinedleastlength;
  extp->logdiagbandwidth = logdiagbandwidth;
  extp->mincoverage = mincoverage;
  extp->display_flag = display_flag;
  extp->use_apos = use_apos;
  extp->xdropbelowscore = xdropbelowscore;
  extp->extendgreedy = extendgreedy;
  extp->extendxdrop = extendxdrop;
  extp->maxalignedlendifference = maxalignedlendifference;
  extp->history_size = history_size;
  extp->perc_mat_history = perc_mat_history;
  extp->extend_char_access = extend_char_access;
  extp->sensitivity = sensitivity;
  extp->matchscore_bias = matchscore_bias;
  extp->weakends = weakends;
  extp->benchmark = benchmark;
  extp->alignmentwidth = alignmentwidth;
  extp->always_polished_ends = always_polished_ends;
  return extp;
}

void gt_diagbandseed_extend_params_delete(GtDiagbandseedExtendParams *extp)
{
  if (extp != NULL) {
    gt_free(extp);
  }
}

/* * * * * K-MER LIST CREATION * * * * */

/* Estimate the number of k-mers in the given encseq. */
static GtUword gt_seed_extend_numofkmers(const GtEncseq *encseq,
                                         unsigned int seedlength,
                                         const GtRange *seqrange)
{
  GtUword lastpos, numofpos, subtract, ratioofspecial;

  const GtUword totalnumofspecial = gt_encseq_specialcharacters(encseq),
                totalnumofpos = gt_encseq_total_length(encseq),
                firstpos = gt_encseq_seqstartpos(encseq, seqrange->start),
                numofseq = gt_range_length(seqrange);
  lastpos = (seqrange->end + 1 == gt_encseq_num_of_sequences(encseq)
             ? totalnumofpos
             : gt_encseq_seqstartpos(encseq, seqrange->end + 1) - 1);
  gt_assert(lastpos >= firstpos);
  numofpos = lastpos - firstpos;

  subtract = MIN(seedlength - 1, gt_encseq_min_seq_length(encseq)) + 1;
  gt_assert(numofpos + 1 >= numofseq * subtract);
  ratioofspecial = MIN(totalnumofspecial * numofpos / totalnumofpos, numofpos);
  return numofpos - MAX(numofseq * subtract - 1, ratioofspecial);
}

/* Returns the position of the next separator following specialrange.start.
   If the end of the encseq is reached, the position behind is returned. */
static GtUword gt_diagbandseed_update_separatorpos(GtRange *specialrange,
                                                   GtSpecialrangeiterator *sri,
                                                   const GtEncseq *encseq,
                                                   GtReadmode readmode)
{
  gt_assert(sri != NULL && specialrange != NULL && encseq != NULL);
  do {
    GtUword idx;
    for (idx = specialrange->start; idx < specialrange->end; idx++) {
      if (gt_encseq_position_is_separator(encseq, idx, readmode)) {
        specialrange->start = idx + 1;
        return idx;
      }
    }
  } while (gt_specialrangeiterator_next(sri, specialrange));
  return gt_encseq_total_length(encseq);
}

/* Add given code and its seqnum and position to a kmer list. */
static void gt_diagbandseed_processkmercode(void *prockmerinfo,
                                            bool firstinrange,
                                            GtUword startpos,
                                            GtCodetype code)
{
  const GtUword array_incr = 256;
  GtDiagbandseedProcKmerInfo *pkinfo;
  GtDiagbandseedKmerPos *kmerposptr = NULL;

  gt_assert(prockmerinfo != NULL);
  pkinfo = (GtDiagbandseedProcKmerInfo *) prockmerinfo;
  GT_GETNEXTFREEINARRAY(kmerposptr,
                        pkinfo->list,
                        GtDiagbandseedKmerPos,
                        array_incr + 0.2 *
                        pkinfo->list->allocatedGtDiagbandseedKmerPos);

  /* check separator positions and determine next seqnum and endpos */
  if (firstinrange) {
    const GtUword endpos = startpos + pkinfo->seedlength - 1;
    while (endpos >= pkinfo->next_separator) {
      pkinfo->seqnum++;
      pkinfo->prev_separator = pkinfo->next_separator + 1;
      pkinfo->next_separator
        = gt_diagbandseed_update_separatorpos(pkinfo->specialrange,
                                              pkinfo->sri,
                                              pkinfo->encseq,
                                              pkinfo->readmode);
      gt_assert(pkinfo->next_separator >= pkinfo->prev_separator);
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
  kmerposptr->code = (pkinfo->readmode == GT_READMODE_FORWARD
                      ? code : gt_kmercode_reverse(code, pkinfo->seedlength));
  /* save endpos and seqnum */
  gt_assert(pkinfo->endpos != UINT_MAX);
  kmerposptr->endpos = pkinfo->endpos;
  pkinfo->endpos = (pkinfo->readmode == GT_READMODE_FORWARD
                    ? pkinfo->endpos + 1 : pkinfo->endpos - 1);
  kmerposptr->seqnum = pkinfo->seqnum;
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
  position = gt_encseq_seqstartpos(pkinfo->encseq, pkinfo->seqnum);
  kc_iter = gt_kmercodeiterator_encseq_new(pkinfo->encseq,
                                           pkinfo->readmode,
                                           pkinfo->seedlength,
                                           position);
  if (pkinfo->seedlength <= pkinfo->totallength) {
    maxpos = pkinfo->totallength + 1 - pkinfo->seedlength;
  }

  /* iterate */
  while (position < maxpos) {
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
    position++;
  }
  gt_kmercodeiterator_delete(kc_iter);
}

/* Return a sorted list of k-mers of given seedlength from specified encseq.
 * Only sequences in seqrange will be taken into account.
 * The caller is responsible for freeing the result. */
GtArrayGtDiagbandseedKmerPos gt_diagbandseed_get_kmers(const GtEncseq *encseq,
                                                       unsigned int seedlength,
                                                       GtReadmode readmode,
                                                       const GtRange *seqrange,
                                                       bool debug_kmer,
                                                       bool verbose,
                                                       GtUword known_size,
                                                       FILE *stream)
{
  GtArrayGtDiagbandseedKmerPos list;
  GtDiagbandseedProcKmerInfo pkinfo;
  GtRadixsortinfo *rdxinfo;
  GtRange specialrange;
  GtTimer *timer = NULL;
  GtUword listlen = known_size;

  gt_assert(encseq != NULL);
  gt_assert(seqrange->start <= seqrange->end);
  gt_assert(seqrange->end < gt_encseq_num_of_sequences(encseq));

  if (known_size == 0) {
    listlen = gt_seed_extend_numofkmers(encseq, seedlength, seqrange);
  }

  if (verbose) {
    timer = gt_timer_new();
    fprintf(stream, "# Start fetching %u-mers (expect " GT_WU ")...\n",
            seedlength, listlen);
    gt_timer_start(timer);
  }

  GT_INITARRAY(&list, GtDiagbandseedKmerPos);
  GT_CHECKARRAYSPACEMULTI(&list, GtDiagbandseedKmerPos, listlen);

  pkinfo.list = &list;
  pkinfo.seqnum = seqrange->start;
  pkinfo.endpos = 0;
  pkinfo.encseq = encseq;
  pkinfo.seedlength = seedlength;
  pkinfo.readmode = readmode;
  if (seqrange->end + 1 == gt_encseq_num_of_sequences(encseq)) {
    pkinfo.totallength = gt_encseq_total_length(encseq);
  } else {
    /* start position of following sequence, minus separator position */
    pkinfo.totallength = gt_encseq_seqstartpos(encseq, seqrange->end + 1) - 1;
  }
  pkinfo.prev_separator = gt_encseq_seqstartpos(encseq, seqrange->start);
  if (gt_encseq_has_specialranges(encseq)) {
    bool search = true;
    pkinfo.sri = gt_specialrangeiterator_new(encseq, true);
    while (search && gt_specialrangeiterator_next(pkinfo.sri, &specialrange)) {
      search = specialrange.end < pkinfo.prev_separator ? true : false;
    }
    specialrange.start = pkinfo.prev_separator;
    pkinfo.specialrange = &specialrange;
    pkinfo.next_separator
      = gt_diagbandseed_update_separatorpos(pkinfo.specialrange,
                                            pkinfo.sri,
                                            pkinfo.encseq,
                                            pkinfo.readmode);
  } else {
    pkinfo.sri = NULL;
    pkinfo.specialrange = NULL;
    pkinfo.next_separator = pkinfo.totallength;
  }

  if (gt_encseq_has_twobitencoding(encseq) && gt_encseq_wildcards(encseq) == 0)
  {
    /* Use fast access to encseq, requires 2bit-enc and absence of wildcards. */
    getencseqkmers_twobitencoding_slice(encseq,
                                        readmode,
                                        seedlength,
                                        seedlength,
                                        false,
                                        gt_diagbandseed_processkmercode,
                                        (void *) &pkinfo,
                                        NULL,
                                        NULL,
                                        pkinfo.prev_separator,
                                        pkinfo.totallength);
  } else {
    /* Use GtKmercodeiterator for encseq access */
    gt_diagbandseed_get_kmers_kciter(&pkinfo);
  }
  if (gt_encseq_has_specialranges(encseq)) {
    gt_specialrangeiterator_delete(pkinfo.sri);
  }
  listlen = list.nextfreeGtDiagbandseedKmerPos;

  /* reduce size of array to number of entries */
  /* list.allocatedGtDiagbandseedKmerPos = listlen;
  gt_realloc(list.spaceGtDiagbandseedKmerPos,
             listlen * sizeof (GtDiagbandseedKmerPos)); */

  if (debug_kmer) {
    GtDiagbandseedKmerPos *idx = list.spaceGtDiagbandseedKmerPos;
    GtDiagbandseedKmerPos *end = idx + listlen;
    while (idx < end) {
      fprintf(stream, "# Kmer (" GT_LX ",%d,%d)\n",
              idx->code, idx->endpos, idx->seqnum);
      idx++;
    }
  }

  if (verbose) {
    fprintf(stream, "# ...found " GT_WU " %u-mers ", listlen, seedlength);
    gt_timer_show_formatted(timer, "in " GT_WD ".%06ld seconds.\n", stream);
    gt_timer_start(timer);
  }

  /* sort list */
  rdxinfo = gt_radixsort_new_ulongpair(listlen);
  gt_radixsort_inplace_GtUwordPair((GtUwordPair *)
                                   list.spaceGtDiagbandseedKmerPos,
                                   listlen);
  gt_radixsort_delete(rdxinfo);

  if (verbose) {
    fprintf(stream, "# ...sorted " GT_WU " %u-mers ", listlen, seedlength);
    gt_timer_show_formatted(timer, "in " GT_WD ".%06ld seconds.\n", stream);
    gt_timer_delete(timer);
  }

  return list;
}

/* * * * * SEEDPAIR LIST CREATION * * * * */

/* Evaluate the results of the seed pair count histogram */
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
        if (histogram != NULL && !selfcomp) {
          histogram[frequency - 1] += (GtUword)(asegm_end - aptr) *
                                      (GtUword)(bsegm_end - bptr);
        } else {
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
                           GtArrayGtDiagbandseedSeedPair *mlist,
                           unsigned int seedlength,
                           bool reverse,
                           bool verbose,
                           FILE *stream,
                           GtError *err) {
  GtTimer *timer = gt_timer_new();
  GtDiagbandseedSeedPair *curr_sp, *max_sp;
  char *buf1 = gt_malloc(3 * (seedlength + 1) * sizeof *buf1);
  char *buf2 = buf1 + 1 + seedlength;
  char *buf3 = buf2 + 1 + seedlength;
  buf1[seedlength] = buf2[seedlength] = buf3[seedlength] = '\0';

  if (verbose) {
    fprintf(stream, "# Start verifying seed pairs...\n");
    gt_timer_start(timer);
  }

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
        gt_timer_delete(timer);
        GT_FREEARRAY(mlist, GtDiagbandseedSeedPair);
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
        gt_timer_delete(timer);
        GT_FREEARRAY(mlist, GtDiagbandseedSeedPair);
        return -1;
      }
    }
    curr_sp++;
  }
  if (verbose) {
    fprintf(stream, "# ...successfully verified each seed pair in ");
    gt_timer_show_formatted(timer, GT_WD ".%06ld seconds.\n", stream);
  }
  gt_free(buf1);
  gt_timer_delete(timer);
  return 0;
}

/* Return estimated length of mlist, and maxfreq w.r.t. the given memlimit */
int gt_diagbandseed_get_mlen_maxfreq(GtUword *mlen,
                                     GtUword *maxfreq,
                                     GtArrayGtDiagbandseedKmerPos *alist,
                                     GtArrayGtDiagbandseedKmerPos *blist,
                                     GtUword memlimit,
                                     unsigned int endposdiff,
                                     bool norev,
                                     bool nofwd,
                                     bool extend_last,
                                     bool selfcomp,
                                     bool alist_blist_id,
                                     bool verbose,
                                     FILE *stream,
                                     GtError *err)
{
  const GtUword maxgram = MIN(*maxfreq, 8190) + 1; /* Cap on k-mer count */
  const GtUword alen = alist->nextfreeGtDiagbandseedKmerPos;
  const GtUword blen = blist->nextfreeGtDiagbandseedKmerPos;
  const bool both_strands = (norev || nofwd) ? false : true;
  GtUword *histogram = NULL;
  GtTimer *timer = NULL;
  GtUword len_used;
  int had_err = 0;

  if (memlimit == GT_UWORD_MAX) {
    return 0; /* no histogram calculation */
  }

  if (verbose) {
    timer = gt_timer_new();
    fprintf(stream, "# Start calculating k-mer frequency histogram...\n");
    gt_timer_start(timer);
  }

  len_used = alen;
  if (!selfcomp || !norev) {
    len_used += blen;
  }
  if (!selfcomp && both_strands && extend_last) {
    len_used += blen;
  }

  /* build histogram; histogram[maxgram] := estimation for mlen */
  histogram = gt_calloc(maxgram + 1, sizeof *histogram);
  gt_diagbandseed_merge(NULL, /* mlist not needed: just count */
                        alist->spaceGtDiagbandseedKmerPos,
                        alen,
                        blist->spaceGtDiagbandseedKmerPos,
                        blen,
                        maxfreq,
                        maxgram,
                        memlimit,
                        histogram,
                        alist_blist_id ? endposdiff : 0,
                        selfcomp,
                        len_used);
  *mlen = histogram[maxgram];
  gt_free(histogram);

  if (verbose) {
    gt_timer_show_formatted(timer,
                            "# ...finished in " GT_WD ".%06ld seconds.\n",
                            stream);
    gt_timer_delete(timer);
  }

  /* check maxfreq value */
  if (*maxfreq == 0 || (*maxfreq == 1 && alist_blist_id)) {
    gt_error_set(err,
                 "option -memlimit too strict: need at least " GT_WU "MB",
                 (*mlen >> 20) + 1);
    *mlen = 0;
    had_err = -1;
    gt_timer_delete(timer);
  } else if (verbose) {
    if (*maxfreq == GT_UWORD_MAX) {
      fprintf(stream, "# Disable k-mer maximum frequency, ");
    } else {
      fprintf(stream, "# Set k-mer maximum frequency to " GT_WU ", ", *maxfreq);
    }
    fprintf(stream, "expect " GT_WU " seed pairs.\n", *mlen);
  } else if (*maxfreq <= 5) {
    gt_warning("Only k-mers occurring <= " GT_WU " times will be considered, "
               "due to small memlimit.", *maxfreq);
  }

  return had_err;
}

/* Return a sorted list of SeedPairs from given alist and blist.
 * Parameter known_size > 0 can be given to allocate the memory beforehand.
 * The caller is responsible for freeing the result. */
GtArrayGtDiagbandseedSeedPair gt_diagbandseed_get_seedpairs(
                                  GtArrayGtDiagbandseedKmerPos *alist,
                                  GtArrayGtDiagbandseedKmerPos *blist,
                                  GtUword maxfreq,
                                  GtUword known_size,
                                  unsigned int endposdiff,
                                  bool selfcomp,
                                  bool delete_blist,
                                  bool debug_seedpair,
                                  bool verbose,
                                  FILE *stream)
{
  GtArrayGtDiagbandseedSeedPair mlist;
  GtTimer *timer = NULL;
  GtUword mlen;

  if (verbose) {
    timer = gt_timer_new();
    if (known_size > 0) {
      fprintf(stream, "# Start building " GT_WU " seed pairs...\n",
                      known_size);
    } else {
      fprintf(stream, "# Start building seed pairs...\n");
    }
    gt_timer_start(timer);
  }

  /* allocate mlist space according to seed pair count */
  GT_INITARRAY(&mlist, GtDiagbandseedSeedPair);
  if (known_size > 0) {
    GT_CHECKARRAYSPACEMULTI(&mlist, GtDiagbandseedSeedPair, known_size);
  }

  /* create mlist */
  gt_diagbandseed_merge(&mlist,
                        alist->spaceGtDiagbandseedKmerPos,
                        alist->nextfreeGtDiagbandseedKmerPos,
                        blist->spaceGtDiagbandseedKmerPos,
                        blist->nextfreeGtDiagbandseedKmerPos,
                        &maxfreq,
                        GT_UWORD_MAX, /* maxgram not needed */
                        GT_UWORD_MAX, /* memlimit not needed */
                        NULL, /* histogram not needed: save seed pairs */
                        endposdiff,
                        selfcomp,
                        0); /* len_used not needed */
  mlen = mlist.nextfreeGtDiagbandseedSeedPair;

  if (delete_blist) {
    GT_FREEARRAY(blist, GtDiagbandseedKmerPos);
    blist = NULL;
  }

  if (verbose) {
    fprintf(stream, "# ...collected " GT_WU " seed pairs ", mlen);
    gt_timer_show_formatted(timer, "in " GT_WD ".%06ld seconds.\n", stream);
  }

  /* sort mlist */
  if (mlen > 0) {
    GtDiagbandseedSeedPair *mspace = mlist.spaceGtDiagbandseedSeedPair;
    GtRadixsortinfo *rdxinfo;

    if (verbose) {
      gt_timer_start(timer);
    }

    rdxinfo = gt_radixsort_new_uint64keypair(mlen);
    gt_radixsort_inplace_Gtuint64keyPair((Gtuint64keyPair*) mspace, mlen);
    gt_radixsort_delete(rdxinfo);

    if (verbose) {
      fprintf(stream, "# ...sorted " GT_WU " seed pairs in ", mlen);
      gt_timer_show_formatted(timer, GT_WD ".%06ld seconds.\n", stream);
    }

    if (debug_seedpair) {
      GtDiagbandseedSeedPair *curr_sp;
      for (curr_sp = mspace; curr_sp < mspace + mlen; curr_sp++) {
        fprintf(stream, "# SeedPair (%d,%d,%d,%d)\n", curr_sp->aseqnum,
                curr_sp->bseqnum, curr_sp->apos, curr_sp->bpos);
      }
    }
  }

  if (verbose) {
    gt_timer_delete(timer);
  }
  return mlist;
}

/* * * * * SEED EXTENSION * * * * */

/* start seed extension for seed pairs in mlist */
static void gt_diagbandseed_process_seeds(GtArrayGtDiagbandseedSeedPair *mlist,
                                          const GtDiagbandseedExtendParams *arg,
                                          void *processinfo,
                                          GtQuerymatchoutoptions *querymoutopt,
                                          const GtEncseq *aencseq,
                                          const GtEncseq *bencseq,
                                          unsigned int seedlength,
                                          bool reverse,
                                          bool verbose,
                                          FILE *stream)
{
  GtDiagbandseedScore *score = NULL;
  GtDiagbandseedPosition *lastp = NULL;
  GtExtendSelfmatchRelativeFunc extend_selfmatch_relative_function = NULL;
  GtExtendQuerymatchRelativeFunc extend_querymatch_relative_function = NULL;
  GtProcessinfo_and_querymatchspaceptr info_querymatch;
  const GtUword amaxlen = gt_encseq_max_seq_length(aencseq),
                bmaxlen = gt_encseq_max_seq_length(bencseq);
  const GtDiagbandseedSeedPair *lm = NULL, *maxsegm, *nextsegm, *idx;
  const GtUword ndiags = (amaxlen >> arg->logdiagbandwidth) +
                         (bmaxlen >> arg->logdiagbandwidth) + 2;
  const GtUword minsegmentlen = (arg->mincoverage - 1) / seedlength + 1;
  GtUword mlen = 0, diag = 0;
  bool firstinrange = true;
  GtUword count_extensions = 0;
  const GtReadmode query_readmode = ((arg->extendgreedy || arg->extendxdrop)
                                     && reverse
                                     ? GT_READMODE_REVCOMPL
                                     : GT_READMODE_FORWARD);
  GtTimer *timer = NULL;
#ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
  GtUword *seedhistogram = NULL;
  GtUword seedcount = 0;
#endif

  gt_assert(mlist != NULL);
  mlen = mlist->nextfreeGtDiagbandseedSeedPair; /* mlist length  */
  lm = mlist->spaceGtDiagbandseedSeedPair;      /* mlist pointer */

  if (mlen < minsegmentlen || mlen == 0) {
    GT_FREEARRAY(mlist, GtDiagbandseedSeedPair);
    return;
  }

  /* select extension method */
  info_querymatch.processinfo = processinfo;
  if (arg->extendgreedy) {
    extend_selfmatch_relative_function = gt_greedy_extend_selfmatch_relative;
    extend_querymatch_relative_function = gt_greedy_extend_querymatch_relative;
  } else if (arg->extendxdrop) {
    extend_selfmatch_relative_function = gt_xdrop_extend_selfmatch_relative;
    extend_querymatch_relative_function = gt_xdrop_extend_querymatch_relative;
  } else { /* no seed extension */
    GT_FREEARRAY(mlist, GtDiagbandseedSeedPair);
    return;
  }

  if (verbose) {
    timer = gt_timer_new();
    if (arg->extendgreedy) {
      fprintf(stream, "# Start greedy seed pair extension...\n");
    } else {
      fprintf(stream, "# Start xdrop seed pair extension...\n");
    }
    fprintf(stream, "# Columns: alen aseq astartpos strand blen bseq bstartpos "
            "score editdist identity\n");
    gt_timer_start(timer);
  }

  info_querymatch.querymatchspaceptr = gt_querymatch_new();
  gt_querymatch_display_set(info_querymatch.querymatchspaceptr,
                            arg->display_flag);
  if (querymoutopt != NULL) {
    gt_querymatch_outoptions_set(info_querymatch.querymatchspaceptr,
                                 querymoutopt);
  }
  gt_querymatch_query_readmode_set(info_querymatch.querymatchspaceptr,
                                   query_readmode);
  gt_querymatch_file_set(info_querymatch.querymatchspaceptr, stream);

  /* score[0] and score[ndiags+1] remain zero for boundary */
  score = gt_calloc(ndiags + 2, sizeof *score);
  lastp = gt_calloc(ndiags, sizeof *lastp);
  maxsegm = lm + mlen - minsegmentlen;
  nextsegm = lm;

#ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
  seedhistogram = (GtUword *)gt_calloc(GT_DIAGBANDSEED_SEEDHISTOGRAM,
                                       sizeof *seedhistogram);
#endif

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
              >> arg->logdiagbandwidth;
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
             >> arg->logdiagbandwidth;
      if ((GtUword)MAX(score[diag + 2], score[diag]) + (GtUword)score[diag + 1]
          >= arg->mincoverage)
      {
        /* relative seed start position in A and B */
        const GtUword bstart = (GtUword) (idx->bpos + 1 - seedlength);
        const GtUword astart = (GtUword) (idx->apos + 1 - seedlength);
#ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
        seedcount++;
#endif

        if (firstinrange ||
            !gt_querymatch_overlap(info_querymatch.querymatchspaceptr,
                                   idx->apos, idx->bpos, arg->use_apos))
        {
          /* extend seed */
          const GtQuerymatch *querymatch = NULL;

          if (aencseq == bencseq) {
            querymatch = extend_selfmatch_relative_function(&info_querymatch,
                                                            aencseq,
                                                            idx->aseqnum,
                                                            astart,
                                                            idx->bseqnum,
                                                            bstart,
                                                            seedlength,
                                                            query_readmode);
          } else {
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
          count_extensions++;
          if (querymatch != NULL) {
            firstinrange = false;
            /* show extension results */
            if (gt_querymatch_check_final(querymatch, arg->errorpercentage,
                                          arg->userdefinedleastlength))
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
             >> arg->logdiagbandwidth;
      score[diag + 1] = 0;
      lastp[diag] = 0;
    }

#ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
    seedhistogram[MIN(GT_DIAGBANDSEED_SEEDHISTOGRAM - 1, seedcount)]++;
    seedcount = 0;
#endif
  }
  gt_querymatch_delete(info_querymatch.querymatchspaceptr);
  gt_free(score);
  gt_free(lastp);

#ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
  fprintf(stream, "# seed histogram:");
  for (seedcount = 0; seedcount < GT_DIAGBANDSEED_SEEDHISTOGRAM; seedcount++) {
    if (seedcount % 10 == 0) {
      fprintf(stream, "\n#\t");
    }
    fprintf(stream, GT_WU "\t", seedhistogram[seedcount]);
  }
  fprintf(stream, "\n");
  gt_free(seedhistogram);
#endif
  if (verbose) {
    fprintf(stream, "# ...finished " GT_WU " seed pair extension%s ",
            count_extensions, count_extensions > 1 ? "s" : "");
    gt_timer_show_formatted(timer, "in " GT_WD ".%06ld seconds.\n", stream);
    gt_timer_delete(timer);
  }
  GT_FREEARRAY(mlist, GtDiagbandseedSeedPair);
}

/* * * * * ALGORITHM STEPS * * * * */

/* Go through the different steps of the seed and extend algorithm. */
static int gt_diagbandseed_algorithm(const GtDiagbandseedInfo *arg,
                                     const GtRange *aseqrange,
                                     const GtRange *bseqrange,
                                     GtArrayGtDiagbandseedKmerPos *alist,
                                     FILE *stream,
                                     GtError *err)
{
  GtArrayGtDiagbandseedKmerPos blist, clist;
  GtArrayGtDiagbandseedSeedPair mlist, mrevlist;
  GtUword alen, blen = 0, mlen = 0, mrevlen = 0, maxfreq;
  unsigned int endposdiff;
  int had_err = 0;
  bool alist_blist_id, both_strands, selfcomp, equalranges;
  GtDiagbandseedExtendParams *extp = NULL;
  Polishing_info *pol_info = NULL;
  void *processinfo = NULL;
  GtQuerymatchoutoptions *querymoutopt = NULL;

  gt_assert(arg != NULL && alist != NULL);
  gt_assert(aseqrange != NULL && bseqrange != NULL);

  maxfreq = arg->maxfreq;
  endposdiff = !arg->overlappingseeds ? arg->seedlength : 1;
  selfcomp = (arg->bencseq == arg->aencseq &&
              gt_range_overlap(aseqrange, bseqrange))
              ? true : false;
  equalranges = gt_range_compare(aseqrange, bseqrange) == 0 ? true : false;
  alist_blist_id = selfcomp && !arg->nofwd && equalranges ? true : false;
  both_strands = (arg->norev || arg->nofwd) ? false : true;
  alen = alist->nextfreeGtDiagbandseedKmerPos;

  /* Second k-mer list */
  if (alist_blist_id) {
    /* compare reads of encseq A with themselves */
    blist = *alist;
    blen = alen;
  } else {
    const GtReadmode readmode = arg->nofwd ? GT_READMODE_COMPL
                                           : GT_READMODE_FORWARD;
    const GtUword known_size = (selfcomp && equalranges) ? alen : 0;
    blist = gt_diagbandseed_get_kmers(arg->bencseq,
                                      arg->seedlength,
                                      readmode,
                                      bseqrange,
                                      arg->debug_kmer,
                                      arg->verbose,
                                      known_size,
                                      stream);
    blen = blist.nextfreeGtDiagbandseedKmerPos;
  }

  had_err = gt_diagbandseed_get_mlen_maxfreq(&mlen,
                                             &maxfreq,
                                             alist,
                                             &blist,
                                             arg->memlimit,
                                             endposdiff,
                                             arg->norev,
                                             arg->nofwd,
                                             arg->extend_last,
                                             selfcomp,
                                             alist_blist_id,
                                             arg->verbose,
                                             stream,
                                             err);
  if (had_err && !alist_blist_id) {
    GT_FREEARRAY(&blist, GtDiagbandseedKmerPos);
  }

  if (!had_err) {
    mlist = gt_diagbandseed_get_seedpairs(alist,
                                          &blist,
                                          maxfreq,
                                          mlen,
                                          alist_blist_id ? endposdiff : 0,
                                          selfcomp,
                                          !alist_blist_id, /* delete_blist */
                                          arg->debug_seedpair,
                                          arg->verbose,
                                          stream);
    mlen = mlist.nextfreeGtDiagbandseedSeedPair;

    if (arg->verify && mlen > 0) {
      had_err = gt_diagbandseed_verify(arg->aencseq,
                                       arg->bencseq,
                                       &mlist,
                                       arg->seedlength,
                                       arg->nofwd,
                                       arg->verbose,
                                       stream,
                                       err);
    }
  }

  if (had_err) {
    return had_err;
  }

  /* Create extension info objects */
  extp = arg->extp;
  if (extp->extendgreedy) {
    GtGreedyextendmatchinfo *grextinfo = NULL;
    const double weak_errorperc = (double)(extp->weakends
                                           ? MAX(extp->errorpercentage, 20)
                                           : extp->errorpercentage);

    pol_info = polishing_info_new_with_bias(weak_errorperc,
                                            extp->matchscore_bias,
                                            extp->history_size);
    grextinfo = gt_greedy_extend_matchinfo_new(extp->errorpercentage,
                                               extp->maxalignedlendifference,
                                               extp->history_size,
                                               extp->perc_mat_history,
                                               extp->userdefinedleastlength,
                                               extp->extend_char_access,
                                               extp->sensitivity,
                                               pol_info);
    if (extp->benchmark) {
      gt_greedy_extend_matchinfo_silent_set(grextinfo);
    }
    processinfo = (void *)grextinfo;
  } else if (extp->extendxdrop) {
    GtXdropmatchinfo *xdropinfo = NULL;
    gt_assert(extp->extendgreedy == false);
    xdropinfo = gt_xdrop_matchinfo_new(extp->userdefinedleastlength,
                                       extp->errorpercentage,
                                       extp->xdropbelowscore,
                                       extp->sensitivity);
    if (extp->benchmark) {
      gt_xdrop_matchinfo_silent_set(xdropinfo);
    }
    processinfo = (void *)xdropinfo;
  }
  if (extp->extendxdrop || extp->alignmentwidth > 0) {
    querymoutopt = gt_querymatchoutoptions_new(true,
                                               false,
                                               extp->alignmentwidth);
    if (extp->extendxdrop || extp->extendgreedy) {
      const GtUword sensitivity = extp->extendxdrop ? 100UL : extp->sensitivity;
      gt_querymatchoutoptions_extend(querymoutopt,
                                     extp->errorpercentage,
                                     extp->maxalignedlendifference,
                                     extp->history_size,
                                     extp->perc_mat_history,
                                     extp->extend_char_access,
                                     extp->weakends,
                                     sensitivity,
                                     extp->matchscore_bias,
                                     extp->always_polished_ends,
                                     extp->display_flag);
    }
  }

  /* process first mlist, unless extend_last */
  if (!arg->extend_last) {
    gt_diagbandseed_process_seeds(&mlist,
                                  arg->extp,
                                  processinfo,
                                  querymoutopt,
                                  arg->aencseq,
                                  arg->bencseq,
                                  arg->seedlength,
                                  arg->nofwd,
                                  arg->verbose,
                                  stream);
  }

  /* Third (reverse) k-mer list */
  if (both_strands) {
    clist = gt_diagbandseed_get_kmers(arg->bencseq,
                                      arg->seedlength,
                                      GT_READMODE_COMPL,
                                      bseqrange,
                                      arg->debug_kmer,
                                      arg->verbose,
                                      blen,
                                      stream);
    gt_assert(blen == clist.nextfreeGtDiagbandseedKmerPos);

    had_err = gt_diagbandseed_get_mlen_maxfreq(&mrevlen,
                                               &maxfreq,
                                               alist,
                                               &clist,
                                               arg->memlimit,
                                               0, /*endposdiff */
                                               arg->norev,
                                               arg->nofwd,
                                               arg->extend_last,
                                               selfcomp,
                                               alist_blist_id,
                                               arg->verbose,
                                               stream,
                                               err);
    if (had_err) {
      GT_FREEARRAY(&clist, GtDiagbandseedKmerPos);
    } else {
      mrevlist = gt_diagbandseed_get_seedpairs(alist,
                                               &clist,
                                               maxfreq,
                                               mrevlen,
                                               0, /* endposdiff */
                                               selfcomp,
                                               true, /* delete_clist */
                                               arg->debug_seedpair,
                                               arg->verbose,
                                               stream);
      mrevlen = mrevlist.nextfreeGtDiagbandseedSeedPair;

      if (arg->verify && mrevlen > 0) {
        had_err = gt_diagbandseed_verify(arg->aencseq,
                                         arg->bencseq,
                                         &mrevlist,
                                         arg->seedlength,
                                         true,
                                         arg->verbose,
                                         stream,
                                         err);
      }
    }
  }

  /* Process first mlist, unless already done */
  if (arg->extend_last) {
    if (!had_err) {
      gt_diagbandseed_process_seeds(&mlist,
                                    arg->extp,
                                    processinfo,
                                    querymoutopt,
                                    arg->aencseq,
                                    arg->bencseq,
                                    arg->seedlength,
                                    arg->nofwd,
                                    arg->verbose,
                                    stream);
    } else {
      GT_FREEARRAY(&mlist, GtDiagbandseedSeedPair);
    }
  }

  /* Process second (reverse) mlist */
  if (!had_err && both_strands) {
    gt_diagbandseed_process_seeds(&mrevlist,
                                  arg->extp,
                                  processinfo,
                                  querymoutopt,
                                  arg->aencseq,
                                  arg->bencseq,
                                  arg->seedlength,
                                  true,
                                  arg->verbose,
                                  stream);
  }

  /* Clean up */
  if (extp->extendgreedy) {
    polishing_info_delete(pol_info);
    gt_greedy_extend_matchinfo_delete((GtGreedyextendmatchinfo *)processinfo);
  } else if (extp->extendxdrop) {
    gt_xdrop_matchinfo_delete((GtXdropmatchinfo *)processinfo);
  }
  if (extp->extendxdrop || extp->alignmentwidth > 0) {
    gt_querymatchoutoptions_delete(querymoutopt);
  }
  return had_err;
}

#ifdef GT_THREADS_ENABLED
typedef struct{
  const GtDiagbandseedInfo *arg;
  GtArrayGtDiagbandseedKmerPos *alist;
  FILE *stream;
  GtError *err;
  const GtRange *aseqrange;
  const GtRange *bseqrange_start;
  const GtRange *bseqrange_end;
  int *had_err;
}GtDiagbandseedThreadInfo;

static void gt_diagbandseed_thread_info_set(GtDiagbandseedThreadInfo *ti,
                                            const GtDiagbandseedInfo *arg,
                                            GtArrayGtDiagbandseedKmerPos *alist,
                                            FILE *stream,
                                            GtError *err,
                                            int *had_err,
                                            const GtRange *aseqrange,
                                            const GtRange *bseqrange_start,
                                            const GtRange *bseqrange_end,
                                            GtUword num_runs_per_thread)
{
  ti->arg = arg;
  ti->alist = alist;
  ti->stream = stream;
  ti->err = err;
  ti->aseqrange = aseqrange;
  ti->bseqrange_start = bseqrange_start;
  ti->bseqrange_end = MIN(bseqrange_start + num_runs_per_thread,
                          bseqrange_end);
  ti->had_err = had_err;
}

static void *gt_diagbandseed_thread_algorithm(void *thread_info)
{
  GtDiagbandseedThreadInfo *info = (GtDiagbandseedThreadInfo *)thread_info;
  const GtRange *bseqrange = info->bseqrange_start;

  while (!(*info->had_err) && bseqrange < info->bseqrange_end) {
    *info->had_err = gt_diagbandseed_algorithm(info->arg,
                                               info->aseqrange,
                                               bseqrange,
                                               info->alist,
                                               info->stream,
                                               info->err);
    bseqrange++;
  }
  return NULL;
}
#endif

/* Run the algorithm by iterating over all combinations of sequence ranges. */
int gt_diagbandseed_run(const GtDiagbandseedInfo *arg,
                        const GtRange *aseqranges,
                        const GtRange *bseqranges,
                        GtUword anumseqranges,
                        GtUword bnumseqranges,
                        GtError *err)
{
  const bool self = arg->aencseq == arg->bencseq ? true : false;
  GtArrayGtDiagbandseedKmerPos alist;
  GtUword aidx, bidx;
  int had_err = 0;
#ifdef GT_THREADS_ENABLED
  GtDiagbandseedThreadInfo *tinfo = gt_malloc(gt_jobs * sizeof *tinfo);
  int *t_errors = gt_calloc(gt_jobs, sizeof *t_errors);
  FILE *stream[gt_jobs];
  GtStrArray *files = gt_str_array_new();
  GtStr *str = gt_str_new();
  unsigned int tidx;

  /* create output streams */
  stream[0] = stdout;
  for (tidx = 1; tidx < gt_jobs; tidx++) {
    gt_str_append_cstr(str, "thread");
    gt_str_append_uint(str, tidx);
    gt_str_append_cstr(str, ".tmp");
    gt_str_array_add(files, str);
    gt_str_reset(str);
    stream[tidx] = gt_xfopen(gt_str_array_get(files, gt_str_array_size(files)
                                              - 1), "w+");
  }
  gt_str_delete(str);
#endif

  for (aidx = 0; !had_err && aidx < anumseqranges; aidx++) {
    /* create alist here to prevent redundant calculations */
    alist = gt_diagbandseed_get_kmers(arg->aencseq,
                                      arg->seedlength,
                                      GT_READMODE_FORWARD,
                                      aseqranges + aidx,
                                      arg->debug_kmer,
                                      arg->verbose,
                                      0,
                                      stdout);
    bidx = self ? aidx : 0;

#ifdef GT_THREADS_ENABLED
    if (gt_jobs <= 1) {
#endif
      while (!had_err && bidx < bnumseqranges) {
        if (arg->verbose && (anumseqranges > 1 || bnumseqranges > 1)) {
          printf("# Process part " GT_WU " vs part " GT_WU "\n",
                 aidx + 1, bidx + 1);
        }
        /* start algorithm with chosen sequence ranges */
        had_err = gt_diagbandseed_algorithm(arg,
                                            aseqranges + aidx,
                                            bseqranges + bidx,
                                            &alist,
                                            stdout,
                                            err);
        bidx++;
      }
#ifdef GT_THREADS_ENABLED
    } else {
      const GtUword num_runs = bnumseqranges - bidx;
      const GtUword num_runs_per_thread = (num_runs - 1) / gt_jobs + 1;
      const GtUword num_threads = (num_runs - 1) / num_runs_per_thread + 1;
      GtArray *threads = gt_array_new(sizeof (GtThread *));
      GtThread *thread;
      int *t_err = t_errors;

      gt_assert(bidx < bnumseqranges);
      gt_assert(num_threads <= gt_jobs);

      /* start additional threads */
      for (tidx = 1; !had_err && tidx < num_threads; tidx++) {
        bidx += num_runs_per_thread;
        gt_diagbandseed_thread_info_set(tinfo + tidx,
                                        arg,
                                        &alist,
                                        stream[tidx],
                                        err,
                                        t_errors + tidx,
                                        aseqranges + aidx,
                                        bseqranges + bidx,
                                        bseqranges + bnumseqranges,
                                        num_runs_per_thread);
        if ((thread = gt_thread_new(gt_diagbandseed_thread_algorithm,
                                    tinfo + tidx, err)) != NULL) {
          gt_array_add(threads, thread);
        } else {
          had_err = -1;
        }
      }

      /* start main thread */
      if (!had_err) {
        bidx = self ? aidx : 0;
        gt_diagbandseed_thread_info_set(tinfo,
                                        arg,
                                        &alist,
                                        stream[0],
                                        err,
                                        t_errors,
                                        aseqranges + aidx,
                                        bseqranges + bidx,
                                        bseqranges + bnumseqranges,
                                        num_runs_per_thread);
        gt_diagbandseed_thread_algorithm(tinfo);
      }

      /* clean up */
      for (tidx = 0; tidx < gt_array_size(threads); tidx++) {
        thread = *(GtThread**) gt_array_get(threads, tidx);
        if (!had_err) gt_thread_join(thread);
        gt_thread_delete(thread);
      }
      gt_array_delete(threads);

      while (!had_err && t_err < t_errors + num_threads) {
        had_err = *(t_err++);
      }
    }
#endif
    GT_FREEARRAY(&alist, GtDiagbandseedKmerPos);
  }
#ifdef GT_THREADS_ENABLED
  gt_free(tinfo);
  gt_free(t_errors);
  for (tidx = 1; tidx < gt_jobs; tidx++) {
    gt_xfclose(stream[tidx]);
  }
  /* print the threads' output to stdout */
  if (!had_err) {
    const size_t block = 512;
    size_t size;
    char buffer[block];
    FILE *stream;
    for (tidx = 0; tidx < gt_str_array_size(files); tidx++) {
      stream = gt_xfopen(gt_str_array_get(files, tidx), "r");
      while ((size = gt_xfread(buffer, sizeof (char), block, stream)) > 0) {
        gt_xfwrite(buffer, sizeof (char), size, stdout);
      }
      gt_xfclose(stream);
    }
  }
  /* remove temporary files */
  for (tidx = 0; tidx < gt_str_array_size(files); tidx++) {
    if (remove(gt_str_array_get(files, tidx)) != 0) {
      /* ignore failure */
    }
  }
  gt_str_array_delete(files);
#endif
  return had_err;
}
