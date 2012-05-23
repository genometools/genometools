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
  const GtEncseq *encseq;
  unsigned int alphasize;
  unsigned long firstmirrorpos;
  unsigned long totallength;

  unsigned int k;
  unsigned int c;
  unsigned int *count;
  unsigned long *kpositions;

  unsigned int currentchar;
  bool seprange;
  bool alltrusted;

  /* stats */
  unsigned long nofkmergroups,
                nofkmeritvs,
                nofkmers,
                nofcorrections;
};

static inline void gt_randomcodes_correct_data_reset(
    GtRandomcodesCorrectData *data)
{
  unsigned int i;
  data->seprange = false;
  data->alltrusted = true;
  data->currentchar = 0;
  for (i = 0; i < data->alphasize; i++)
    data->count[i] = 0;
}

static inline void gt_randomcodes_correct_process_kmer_itv(
    const unsigned long *suffixes,
    unsigned long nofsuffixes,
    GtRandomcodesCorrectData *cdata)
{
  unsigned int i, c = cdata->c;
  cdata->nofkmeritvs++;
  cdata->nofkmers += nofsuffixes;
  if (nofsuffixes < cdata->c)
  {
    cdata->alltrusted = false;
    c = nofsuffixes;
  }
  for (i = 0; i < c; i++)
  {
    (cdata->kpositions + cdata->currentchar * cdata->c)[i] = suffixes[i];
  }
  cdata->count[cdata->currentchar] = (unsigned int)nofsuffixes;
  cdata->currentchar++;
}

static inline void gt_randomcodes_correct_process_kmergroup_end(
    const GtSeqnumrelpos *snrp, GtRandomcodesCorrectData *cdata)
{
  cdata->nofkmergroups++;
  if (!cdata->alltrusted)
  {
    unsigned int i, max_count = 0, countpos = 0;
    for (i = 0; i < cdata->alphasize; i++)
    {
      if (cdata->count[i] > max_count)
      {
        max_count = cdata->count[i];
        countpos = i;
      }
    }
    if (max_count >= cdata->c)
    {
      const unsigned long trusted_char_seqnumrelpos =
        (cdata->kpositions + countpos * cdata->c)[0];
      GtUchar trusted_char = gt_encseq_get_encoded_char_nospecial(cdata->encseq,
          gt_encseq_seqstartpos(cdata->encseq, gt_seqnumrelpos_decode_seqnum(
              snrp,trusted_char_seqnumrelpos)) + gt_seqnumrelpos_decode_relpos(
              snrp,trusted_char_seqnumrelpos) + cdata->k - 1,
          GT_READMODE_FORWARD);
      for (i = 0; i < cdata->alphasize; i++)
      {
        unsigned int j;
        if (cdata->count[i] < cdata->c)
        {
          for (j = 0; j < cdata->count[i]; j++)
          {
            const unsigned long
              seqnumrelpos = (cdata->kpositions + i * cdata->c)[j];
            unsigned long abspos;
            abspos = gt_encseq_seqstartpos(cdata->encseq,
                gt_seqnumrelpos_decode_seqnum(snrp,seqnumrelpos)) +
              gt_seqnumrelpos_decode_relpos(snrp,seqnumrelpos) + cdata->k - 1;
            GtUchar newchar = trusted_char;
            if (abspos >= cdata->firstmirrorpos)
            {
              abspos = cdata->totallength - 1UL - abspos;
              newchar = (GtUchar)3 - newchar;
            }
            /*gt_log_log("correct: position %lu char %u", abspos, newchar);*/
            cdata->nofcorrections++;
          }
        }
      }
    }
  }
  gt_randomcodes_correct_data_reset(cdata);
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
          next_itvstart - itvstart, data);
      itvstart = next_itvstart;
      if (lcpvalue < correction_kmersize - 1)
      {
        gt_randomcodes_correct_process_kmergroup_end(snrp, data);
      }
    }
  }
  gt_randomcodes_correct_process_kmer_itv(bucketofsuffixes + itvstart,
      numberofsuffixes - itvstart, data);
  gt_randomcodes_correct_process_kmergroup_end(snrp, data);
  return haserr ? -1 : 0;
}

GtRandomcodesCorrectData *gt_randomcodes_correct_data_new(GtEncseq *encseq,
    unsigned int k, unsigned int c)
{
  GtRandomcodesCorrectData *cdata = gt_malloc(sizeof *cdata);
  cdata->k = k;
  cdata->c = c;
  cdata->encseq = encseq;
  cdata->alphasize = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  cdata->totallength = gt_encseq_total_length(encseq);
  cdata->firstmirrorpos = cdata->totallength;
  if (gt_encseq_is_mirrored(encseq))
    cdata->firstmirrorpos >>= 1;
  cdata->kpositions = gt_malloc(sizeof (*cdata->kpositions) *
      cdata->alphasize * c);
  cdata->count = gt_malloc(sizeof (*cdata->count) * cdata->alphasize);
  gt_randomcodes_correct_data_reset(cdata);
  /* stats */
  cdata->nofkmergroups = 0;
  cdata->nofkmeritvs = 0;
  cdata->nofkmers = 0;
  cdata->nofcorrections = 0;
  return cdata;
}

void gt_randomcodes_correct_data_delete(
    GtRandomcodesCorrectData *cdata)
{
  gt_free(cdata->kpositions);
  gt_free(cdata->count);
  /* output stats */
  gt_log_log("nofkmergroups = %lu", cdata->nofkmergroups);
  gt_log_log("nofkmeritvs = %lu", cdata->nofkmeritvs);
  gt_log_log("nofkmers = %lu", cdata->nofkmers);
  gt_log_log("nofcorrections = %lu", cdata->nofcorrections);
  gt_free(cdata);
}
