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

#include "core/log_api.h"
#include "match/randomcodes-correct.h"

struct GtRandomcodesCorrectData {
  unsigned long nofkmergroups,
                nofkmeritvs,
                nofkmers;
};

static inline void gt_randomcodes_correct_process_kmer_itv(
    GT_UNUSED const unsigned long *suffixes,
    unsigned long nofsuffixes,
    GT_UNUSED const GtSeqnumrelpos *snrp,
    void *data)
{
  GtRandomcodesCorrectData *cdata = data;
  cdata->nofkmeritvs++;
  cdata->nofkmers += nofsuffixes;
}

static inline void gt_randomcodes_correct_process_kmergroup_end(void *data)
{
  GtRandomcodesCorrectData *cdata = data;
  cdata->nofkmergroups++;
}

int gt_randomcodes_correct_process_bucket(void *data,
    const unsigned long *bucketofsuffixes, const GtSeqnumrelpos *snrp,
    const uint16_t *lcptab_bucket, unsigned long numberofsuffixes,
    unsigned int correction_kmersize, GT_UNUSED GtError *err)
{
  unsigned long itvstart, next_itvstart;
  unsigned int lcpvalue;
  bool haserr = false, firstedgefromroot;

  firstedgefromroot = true;
  for (itvstart = 0, next_itvstart = 1UL; next_itvstart < numberofsuffixes;
      next_itvstart++)
  {
    lcpvalue = (unsigned int) lcptab_bucket[next_itvstart];
    if (lcpvalue < correction_kmersize)
    {
      gt_randomcodes_correct_process_kmer_itv(bucketofsuffixes + itvstart,
          next_itvstart - itvstart, snrp, data);
      itvstart = next_itvstart;
      if (lcpvalue < correction_kmersize - 1)
      {
        gt_randomcodes_correct_process_kmergroup_end(data);
      }
    }
  }
  gt_randomcodes_correct_process_kmer_itv(bucketofsuffixes + itvstart,
      numberofsuffixes - itvstart, snrp, data);
  gt_randomcodes_correct_process_kmergroup_end(data);
  return haserr ? -1 : 0;
}

GtRandomcodesCorrectData *gt_randomcodes_correct_data_new()
{
  GtRandomcodesCorrectData *cdata = gt_malloc(sizeof *cdata);
  cdata->nofkmergroups = 0;
  cdata->nofkmeritvs = 0;
  return cdata;
}

void gt_randomcodes_correct_data_delete(
    GtRandomcodesCorrectData *cdata)
{
  gt_log_log("nofkmergroups = %lu", cdata->nofkmergroups);
  gt_log_log("nofkmeritvs = %lu", cdata->nofkmeritvs);
  gt_log_log("nofkmers = %lu", cdata->nofkmers);
  gt_free(cdata);
}
