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
#include "match/firstcodes-spacelog.h"
#include "match/kmercodes.h"
#include "match/randomsamples.h"

struct GtRandomSamples {
  GtEncseq *encseq;
  unsigned long nofsequences;
  unsigned long totallength;
  unsigned int keysize;
  unsigned long nofsamples;
  unsigned long *samples;
  GtTimer *timer;
  GtFirstcodesspacelog *spacelog;
};

GtRandomSamples* gt_randomsamples_new(GtEncseq *encseq, unsigned int keysize,
    GtTimer *timer)
{
  GtRandomSamples *rs = gt_malloc(sizeof (*rs));
  rs->encseq = encseq;
  rs->nofsequences = gt_encseq_num_of_sequences(encseq);
  rs->totallength = gt_encseq_total_length(encseq);
  gt_assert(keysize > 0);
  gt_assert(keysize <= 32);
  rs->keysize = keysize;
  rs->samples = NULL;
  rs->timer = timer;
  rs->spacelog = gt_firstcodes_spacelog_new();
  GT_FCI_ADDWORKSPACE(rs->spacelog,"GtRandomSamples", sizeof (*rs));
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
    GT_FCI_SUBTRACTADDWORKSPACE(rs->spacelog,"samples", sizeof (*rs->samples) *
        rs->nofsamples);
  }
  gt_log_log("number of unique samples = %lu", storeindex);
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
    gt_log_log("range of positions = [%lu, %lu]", rs->samples[0],
        rs->samples[rs->nofsamples - 1]);
  }
}

void gt_randomsamples_sample(GtRandomSamples *rs,
    unsigned long sampling_factor)
{
  unsigned long i;
  gt_assert(rs != NULL);
  rs->nofsamples = gt_randomsamples_calculate_nofsamples(rs->encseq,
      rs->nofsequences, rs->totallength, rs->keysize, sampling_factor);
  gt_log_log("nofsamples = %lu", rs->nofsamples);
  rs->samples = gt_malloc(sizeof (*rs->samples) * rs->nofsamples);
  GT_FCI_ADDWORKSPACE(rs->spacelog,"samples", sizeof (*rs->samples) *
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
}

void gt_randomsamples_delete(GtRandomSamples *rs)
{
  if (rs->timer != NULL)
    gt_timer_show_progress(rs->timer, "for cleaning up", stdout);
  if (rs->samples != NULL)
  {
    gt_free(rs->samples);
    GT_FCI_SUBTRACTWORKSPACE(rs->spacelog,"samples");
  }
  GT_FCI_SUBTRACTWORKSPACE(rs->spacelog,"GtRandomSamples");
  gt_firstcodes_spacelog_delete(rs->spacelog);
  gt_free(rs);
}
