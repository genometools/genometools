/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/fa.h"
#include "core/log_api.h"
#include "core/xansi_api.h"
#include "match/randomcodes-find-seldom.h"

struct GtRandomcodesFindSeldomData {
  const GtEncseq *encseq;
  GtUword nofreads;
  GtUword mirror_nofseqs;
  GtBitsequence *seldom_reads;
  GtUword k, c;
  GtUword nofseldomkmers;
};

static inline void gt_randomcodes_find_seldom_process_kmer_itv(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesFindSeldomData *sdata)
{
  GtUword i, relpos, seqnum;
  if (nofsuffixes < sdata->c)
  {
    for (i = 0; i < nofsuffixes; i++)
    {
      relpos = gt_seqnumrelpos_decode_relpos(snrp, suffixes[i]);
      seqnum = gt_seqnumrelpos_decode_seqnum(snrp, suffixes[i]);
      if (seqnum >= sdata->nofreads)
        seqnum = sdata->mirror_nofseqs - seqnum - 1UL;
      if (gt_encseq_seqlength(sdata->encseq, seqnum) - relpos >= sdata->k)
      {
        GT_SETIBIT(sdata->seldom_reads, seqnum);
        sdata->nofseldomkmers++;
        /*
        gt_log_log("seldom "GT_WU"-mer, count="GT_WU", relpos="GT_WU
                   ",seqnum="GT_WU"", sdata->k, nofsuffixes, relpos,
                   seqnum);
         */
      }
    }
  }
}

int gt_randomcodes_find_seldom_process_bucket(void *data,
    const GtUword *bucketofsuffixes, const GtSeqnumrelpos *snrp,
    const uint16_t *lcptab_bucket, GtUword numberofsuffixes,
    GT_UNUSED unsigned int sortingdepth, GT_UNUSED GtError *err)
{
  GtUword itvstart, next_itvstart;
  unsigned int lcpvalue;
  GtRandomcodesFindSeldomData *sdata = data;

  for (itvstart = 0, next_itvstart = 1UL; next_itvstart < numberofsuffixes;
      next_itvstart++)
  {
    lcpvalue = (unsigned int) lcptab_bucket[next_itvstart];
    if (lcpvalue < (unsigned int)sdata->k)
    {
      gt_randomcodes_find_seldom_process_kmer_itv(snrp,
          bucketofsuffixes + itvstart, next_itvstart - itvstart, data);
      itvstart = next_itvstart;
    }
  }
  gt_randomcodes_find_seldom_process_kmer_itv(snrp,
      bucketofsuffixes + itvstart, next_itvstart - itvstart, data);
  return 0;
}

GtRandomcodesFindSeldomData *gt_randomcodes_find_seldom_data_new(
    GtEncseq *encseq, unsigned int k, unsigned int c,
    GtBitsequence *seldom_reads)
{
  GtRandomcodesFindSeldomData *sdata = gt_malloc(sizeof *sdata);
  sdata->k = (GtUword)k;
  sdata->c = (GtUword)c;
  sdata->encseq = encseq;
  sdata->mirror_nofseqs = gt_encseq_num_of_sequences(encseq);
  sdata->nofreads = sdata->mirror_nofseqs;
  sdata->seldom_reads = seldom_reads;
  if (gt_encseq_is_mirrored(sdata->encseq))
    sdata->nofreads >>= 1;
  sdata->nofseldomkmers = 0;
  return sdata;
}

void gt_randomcodes_find_seldom_data_collect_stats(
    GtRandomcodesFindSeldomData *sdata, unsigned int threadnum,
    GtUword *nofseldomkmers)
{
  gt_log_log("thread %u: "GT_WU"", threadnum, sdata->nofseldomkmers);
  *nofseldomkmers += sdata->nofseldomkmers;
}

void gt_randomcodes_find_seldom_data_delete(GtRandomcodesFindSeldomData *sdata)
{
  gt_free(sdata);
}
