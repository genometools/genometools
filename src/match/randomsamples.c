/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/encseq.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/radix_sort.h"
#include "core/unused_api.h"
#include "core/showtime.h"
#include "core/spacecalc.h"
#include "core/thread.h"
#include "match/firstcodes-accum.h"
#include "match/firstcodes-buf.h"
#include "match/firstcodes-cache.h"
#include "match/firstcodes-spacelog.h"
#include "match/kmercodes.h"
#include "match/spmsuftab.h"
#include "match/sfx-maprange.h"
#include "match/randomsamples.h"

struct GtRandomSamples {
  GtEncseq *encseq;
  unsigned long nofsequences,
                totallength,
                maxseqlen;
  unsigned int keysize;
  unsigned int correction_kmersize;
  unsigned long nofsamples;
  unsigned long *samples;
  GtFirstcodesspacelog *spacelog;
  GtCodeposbuffer buf;
  GtSpmsuftab *spmsuftab;
  GtError *err;
  unsigned long codebuffer_total;
  unsigned int flushcount;
  GtRadixsortinfo *radixsort_code;
  GtArrayGtIndexwithcode *binsearchcache;
  GtLogger *default_logger, *verbose_logger;
  GtTimer *timer;
};

GtRandomSamples* gt_randomsamples_new(GtEncseq *encseq, GtTimer *timer,
    unsigned int correction_kmersize, GtLogger *default_logger,
    GtLogger *verbose_logger, GtError *err)
{
  unsigned long maxrelpos;
  unsigned int bitsforrelpos, bitsforseqnum;
  GtRandomSamples *rs = gt_malloc(sizeof (*rs));
  gt_error_check(err);
  rs->err = err;
  rs->encseq = encseq;
  rs->nofsequences = gt_encseq_num_of_sequences(encseq);
  rs->totallength = gt_encseq_total_length(encseq);
  rs->maxseqlen = gt_encseq_max_seq_length(encseq);
  rs->samples = NULL;
  rs->timer = timer;
  rs->spacelog = gt_firstcodes_spacelog_new();
  gt_assert(correction_kmersize > 1U);
  rs->correction_kmersize = correction_kmersize;
  gt_log_log("kmersize for correction = %u", rs->correction_kmersize);
  rs->keysize = MIN(rs->correction_kmersize - 1, GT_UNITSIN2BITENC);
  gt_log_log("kmersize for bucket keys = %u", rs->keysize);
  rs->buf.accum_all = true;
  rs->buf.markprefix = NULL;
  rs->buf.marksuffix = NULL;
  rs->buf.nextfree = 0;
  maxrelpos = (rs->maxseqlen > (unsigned long) rs->correction_kmersize)
    ? rs->maxseqlen - (unsigned long) rs->correction_kmersize : 0;
  bitsforrelpos = gt_determinebitspervalue(maxrelpos);
  rs->buf.snrp = gt_seqnumrelpos_new(bitsforrelpos,encseq);
  bitsforseqnum = gt_determinebitspervalue(rs->nofsequences - 1UL);
  if (bitsforseqnum + bitsforrelpos > (unsigned int) GT_INTWORDSIZE)
  {
    gt_seqnumrelpos_delete(rs->buf.snrp);
    rs->buf.snrp = NULL;
    gt_error_set(rs->err,"cannot process encoded sequences with %lu sequences "
                     "of length up to %lu (%u+%u bits)",
                     rs->nofsequences, rs->maxseqlen, bitsforseqnum,
                     bitsforrelpos);
    gt_randomsamples_delete(rs);
    rs = NULL;
  }
  rs->buf.spaceGtUlongPair = NULL;
  rs->buf.spaceGtUlong = NULL;
  GT_FCI_ADDWORKSPACE(rs->spacelog,"encseq",(size_t)
      gt_encseq_sizeofrep(encseq));
  rs->flushcount = 0;
  rs->codebuffer_total = 0;
  rs->default_logger = default_logger;
  rs->verbose_logger = verbose_logger;
  return rs;
}

#define GT_RANDOMSAMPLES_NOFSAMPLES_MIN 2UL

static inline unsigned long gt_randomsamples_calculate_nofsamples(
    const GtEncseq *encseq, unsigned long nofsequences,
    unsigned long totallength, unsigned int keysize,
    unsigned long sampling_factor)
{
  unsigned long nofkmers = totallength,
                nofnonkmers = (keysize + 1) * nofsequences,
                nofsamples;
  gt_log_log("totallength = %lu", nofkmers);
  gt_log_log("nofsequences = %lu", nofsequences);
  if (nofnonkmers < nofkmers)
  {
    nofkmers -= nofnonkmers;
    gt_log_log("nofkmers = %lu", nofkmers);
  }
  else
  {
    gt_assert(gt_encseq_min_seq_length(encseq) <= keysize);
  }
  nofsamples = nofkmers / sampling_factor;
  return MAX(GT_RANDOMSAMPLES_NOFSAMPLES_MIN, nofsamples);
}

static void gt_randomsamples_rmdup(GtRandomSamples *rs)
{
  unsigned long storeindex, readindex, prev;
  prev = rs->samples[0];
  for (readindex = 1UL, storeindex = 1UL; readindex < rs->nofsamples;
      readindex++)
  {
    if (rs->samples[readindex] != prev)
    {
      prev = rs->samples[readindex];
      if (readindex != storeindex)
        rs->samples[storeindex] = prev;
      storeindex++;
    }
  }
  if (storeindex != readindex)
  {
    rs->nofsamples = storeindex;
    rs->samples = gt_realloc(rs->samples, sizeof (*rs->samples) *
        rs->nofsamples);
    GT_FCI_SUBTRACTADDSPLITSPACE(rs->spacelog,"samples", sizeof (*rs->samples) *
        rs->nofsamples);
  }
  gt_log_log("number of unique samples = %lu (%.3f%%)", storeindex,
      (float)storeindex / (float)readindex * 100.0);
}

static void gt_randomsamples_generate_positions(GtRandomSamples *rs,
    unsigned long sampling_factor, bool sort)
{
  unsigned long i, randmax = sampling_factor;
  bool sorted = true;
  randmax = GT_MULT2(sampling_factor) - GT_DIV32(sampling_factor);
  gt_assert(randmax < rs->totallength);
  rs->samples[0] = gt_rand_max(randmax);
  for (i = 1UL; i < rs->nofsamples; i++)
  {
    unsigned long sn, sp, sl, rp;
    rs->samples[i] = rs->samples[i-1];
    while (true)
    {
      rs->samples[i] += gt_rand_max(randmax);
      if (rs->samples[i] >= rs->totallength)
      {
        rs->samples[i] = 0;
        sorted = false;
      }
      if (!gt_encseq_position_is_separator(rs->encseq, rs->samples[i],
            GT_READMODE_FORWARD))
      {
        sn = gt_encseq_seqnum(rs->encseq, rs->samples[i]);
        sp = gt_encseq_seqstartpos(rs->encseq, sn);
        sl = gt_encseq_seqlength(rs->encseq, sn);
        rp = rs->samples[i] - sp;
        if (rp < sl - rs->keysize)
          break;
      }
    }
  }
  if (sort && !sorted)
  {
    if (rs->timer != NULL)
      gt_timer_show_progress(rs->timer, "to sort sampling positions", stdout);
    gt_radixsort_inplace_ulong(rs->samples, rs->nofsamples);
  }
  if (sorted)
  {
    gt_log_log("range of sampling positions = [%lu, %lu]", rs->samples[0],
        rs->samples[rs->nofsamples - 1]);
  }
}

static int gt_randomsamples_allocspace(GtRandomSamples *rs,
                                    unsigned int numofparts,
                                    unsigned long maximumspace,
                                    unsigned long phase2extra,
                                    GtError *err)
{
  if (maximumspace > 0)
  {
    if ((unsigned long) gt_firstcodes_spacelog_total(rs->spacelog) +
                        phase2extra >= maximumspace)
    {
      gt_error_set(err,"already used %.2f MB of memory and need %.2f MB later "
                       "=> cannot compute index in at most %.2f MB",
                       GT_MEGABYTES(gt_firstcodes_spacelog_total(rs->spacelog)),
                       GT_MEGABYTES(phase2extra),
                       GT_MEGABYTES(maximumspace));
      return -1;
    } else
    {
      size_t remainspace = (size_t) maximumspace -
                           (gt_firstcodes_spacelog_total(rs->spacelog) +
                            phase2extra);

      rs->buf.allocated
        = gt_radixsort_max_num_of_entries_ulong(remainspace);
      if (rs->buf.allocated < rs->nofsamples / 16UL)
      {
        rs->buf.allocated = rs->nofsamples / 16UL;
      }
    }
  } else
  {
    if (numofparts == 0)
    {
      rs->buf.allocated = gt_radixsort_max_num_of_entries_ulong(
                             gt_firstcodes_spacelog_total(rs->spacelog) / 5UL);
    } else
    {
      rs->buf.allocated = rs->nofsamples / 5UL;
    }
  }
  if (rs->buf.allocated < 16UL)
  {
    rs->buf.allocated = 16UL;
  }
  return 0;
}

static void gt_randomsamples_accumulatecounts_flush(void *data)
{
  GtRandomSamples *rs = (GtRandomSamples *) data;

  if (rs->buf.nextfree > 0)
  {
    unsigned long foundindex = ULONG_MAX, foundcode;

    gt_assert(rs->samples != NULL);
    rs->codebuffer_total += rs->buf.nextfree;
    gt_radixsort_inplace_sort(rs->radixsort_code, rs->buf.nextfree);
    foundindex = gt_randomsamples_find_accu(&foundcode,
                                         rs->samples,
                                         rs->nofsamples,
                                         rs->binsearchcache,
                                         rs->buf.spaceGtUlong[0]);
    gt_assert(foundindex != ULONG_MAX);
    /*gt_randomcodes_accumulatecounts_merge(&rs->tab,
          rs->samples,
          rs->differentcodes,
          rs->buf.spaceGtUlong,
          rs->buf.spaceGtUlong
          + rs->buf.nextfree - 1,
          foundindex,
          foundcode);*/
    rs->flushcount++;
    rs->buf.nextfree = 0;
  }
}

static int gt_randomsamples_accumulatecounts_run(GtRandomSamples *rs)
{
  bool haserr = false;
  gt_assert(rs != NULL);
  if (rs->timer != NULL)
  {
    gt_timer_show_progress(rs->timer, "to accumulate counts",stdout);
  }
  gt_assert(rs->buf.allocated > 0);
  rs->radixsort_code = gt_radixsort_new_ulong(rs->buf.allocated);
  rs->buf.spaceGtUlong = gt_radixsort_space_ulong(rs->radixsort_code);
  GT_FCI_ADDWORKSPACE(rs->spacelog, "radixsort_code",
                      gt_radixsort_size(rs->radixsort_code));
  rs->buf.fciptr = rs; /* as we need to give rs to the flush function */
  rs->buf.flush_function = gt_randomsamples_accumulatecounts_flush;
  gt_logger_log(rs->verbose_logger,
      "maximum space for accumulating counts %.2f MB",
      GT_MEGABYTES(gt_firstcodes_spacelog_total(rs->spacelog)));
  gt_firstcodes_accum_runkmerscan(rs->encseq, rs->keysize, rs->keysize,
      &rs->buf);
  gt_randomsamples_accumulatecounts_flush(rs);
  gt_log_log("codebuffer_total = %lu (%.3f%% of all suffixes)",
      rs->codebuffer_total,
      100.0 * (double) rs->codebuffer_total / rs->totallength);
  /*if (rs->firstcodehits > 0)
  {
    gt_assert(rs->flushcount > 0);
    gt_logger_log(logger,"firstcodehits=%lu (%.3f%% of all suffixes), "
                         "%u rounds (avg length %lu)",
                         rs->firstcodehits,
                         100.0 * (double) rs->firstcodehits/
                                          gt_encseq_total_length(encseq),
                         rs->flushcount,
                         rs->codebuffer_total/rs->flushcount);
  }*/
  gt_radixsort_delete(rs->radixsort_code);
  rs->radixsort_code = NULL;
  GT_FCI_SUBTRACTWORKSPACE(rs->spacelog, "radixsort_code");
  return haserr ? -1 : 0;
}

int gt_randomsamples_run(GtRandomSamples *rs, unsigned long sampling_factor,
    unsigned int numofparts, unsigned long maximumspace)
{
  bool haserr = false;
  unsigned long i;
  GtSfxmappedrangelist *sfxmrlist = NULL;
#ifdef GT_THREADS_ENABLED
  GT_UNUSED const unsigned int threads = gt_jobs;
#else
  GT_UNUSED const unsigned int threads = 1U;
#endif
  gt_assert(rs != NULL);

  if (rs->maxseqlen < (unsigned long)rs->correction_kmersize)
    return 0;

  sfxmrlist = gt_Sfxmappedrangelist_new();
  rs->nofsamples = gt_randomsamples_calculate_nofsamples(rs->encseq,
      rs->nofsequences, rs->totallength, rs->keysize, sampling_factor);
  gt_log_log("nofsamples = %lu", rs->nofsamples);
  rs->samples = gt_malloc(sizeof (*rs->samples) * rs->nofsamples);
  GT_FCI_ADDSPLITSPACE(rs->spacelog,"samples", sizeof (*rs->samples) *
      rs->nofsamples);

  if (rs->timer != NULL)
    gt_timer_show_progress(rs->timer, "to generate sampling positions", stdout);
  gt_randomsamples_generate_positions(rs, sampling_factor, true);

  if (rs->timer != NULL)
    gt_timer_show_progress(rs->timer, "to collect sample codes", stdout);
  for (i = 0; i < rs->nofsamples; i++)
  {
    const GtTwobitencoding *twobitenc =
      gt_encseq_twobitencoding_export(rs->encseq);
    rs->samples[i] = gt_kmercode_at_position(twobitenc, rs->samples[i],
        rs->keysize);
  }

  if (rs->timer != NULL)
    gt_timer_show_progress(rs->timer, "to sort sample codes", stdout);
  gt_radixsort_inplace_ulong(rs->samples, rs->nofsamples);

  if (rs->timer != NULL)
    gt_timer_show_progress(rs->timer, "to remove duplicated sample codes",
        stdout);
  gt_randomsamples_rmdup(rs);
  if (!haserr)
  {
    if (gt_randomsamples_allocspace(rs, numofparts, maximumspace, 0,
          rs->err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    gt_randomsamples_accumulatecounts_run(rs);
  }
  gt_Sfxmappedrangelist_delete(sfxmrlist);
  return haserr ? -1 : 0;
}

void gt_randomsamples_delete(GtRandomSamples *rs)
{
  if (rs->timer != NULL)
    gt_timer_show_progress(rs->timer, "for cleaning up", stdout);
  if (rs->samples != NULL)
  {
    gt_free(rs->samples);
    GT_FCI_SUBTRACTSPLITSPACE(rs->spacelog,"samples");
  }
  GT_FCI_SUBTRACTWORKSPACE(rs->spacelog,"encseq");
  gt_firstcodes_spacelog_delete(rs->spacelog);
  gt_seqnumrelpos_delete(rs->buf.snrp);
  gt_free(rs);
}
