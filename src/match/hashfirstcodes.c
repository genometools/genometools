/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include "core/codetype.h"
#include "core/encseq.h"
#include "core/showtime.h"
#include "core/unused_api.h"
#include "extended/uint64hashtable.h"
#include "hashfirstcodes.h"
#include "sfx-suffixer.h"

typedef struct
{
  GtUint64hashtable *table;
  unsigned long differentcodes;
  unsigned int kmersize;
  unsigned long *suftab;
  unsigned long finalpsum;
} GtHashfirstcodes;

static void gt_hashfirstcodes(void *processinfo,
                              GT_UNUSED bool firstinrange,
                              GT_UNUSED unsigned long pos,
                              GtCodetype code)
{
  GtHashfirstcodes *hashfirstcodes = (GtHashfirstcodes *) processinfo;

  if (!gt_uint64hashtable_search(hashfirstcodes->table,
                                 (uint64_t) code,
                                 true))
  {
    hashfirstcodes->differentcodes++;
  }
}

static void gt_hashremainingcodes(void *processinfo,
                                  bool firstinrange,
                                  GT_UNUSED unsigned long pos,
                                  GtCodetype code)
{
  GtUint64hashtable *table = (GtUint64hashtable *) processinfo;

  if (!firstinrange)
  {
    (void) gt_uint64hashtable_search(table,(uint64_t) code,false);
  }
}

static void gt_insertallcodeswithhashtable(void *processinfo,
                                           GT_UNUSED bool firstinrange,
                                           unsigned long pos,
                                           GtCodetype code)
{
  GtHashfirstcodes *hashfirstcodes = (GtHashfirstcodes *) processinfo;

  unsigned long idx = gt_uint64hashtable_insertionindex(hashfirstcodes->table,
                                                        (uint64_t) code);
  if (idx != ULONG_MAX)
  {
    gt_assert(idx < hashfirstcodes->finalpsum);
    hashfirstcodes->suftab[idx] = pos;
  }
}

void hashfirstcodes_getencseqkmers_twobitencoding(const GtEncseq *encseq,
                                                  unsigned int kmersize)
{
  GtTimer *timer = NULL;
  GtHashfirstcodes hashfirstcodes;
  unsigned long numofsequences, psum;
  GT_UNUSED unsigned long totallength;
  /*GT_UNUSED unsigned long countsum;*/
  const unsigned int spmopt = 45U;

  if (gt_showtime_enabled())
  {
    timer = gt_timer_new_with_progress_description("insert first codes into "
                                                   "hashtable");
    gt_timer_start(timer);
  }
  numofsequences = gt_encseq_num_of_sequences(encseq);
  totallength = gt_encseq_total_length(encseq);
  hashfirstcodes.table = gt_uint64hashtable_new((size_t) numofsequences);
  hashfirstcodes.differentcodes = 0;
  getencseqkmers_twobitencoding(encseq,
                                GT_READMODE_FORWARD,
                                kmersize,
                                kmersize,
                                true,
                                gt_hashfirstcodes,
                                &hashfirstcodes,
                                NULL,
                                NULL);
  printf("# number of different codes=%lu (%.4f) in %lu sequences\n",
          hashfirstcodes.differentcodes,
          (double) hashfirstcodes.differentcodes/numofsequences,
          numofsequences);
  /*countsum = gt_uint64hashtable_countsum_get(hashfirstcodes.table);*/
  gt_assert(gt_uint64hashtable_countsum_get(hashfirstcodes.table)
               == numofsequences);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "accumulate counts",stdout);
  }
  getencseqkmers_twobitencoding(encseq,
                                GT_READMODE_FORWARD,
                                kmersize,
                                spmopt,
                                false,
                                gt_hashremainingcodes,
                                hashfirstcodes.table,
                                NULL,
                                NULL);
  psum = gt_uint64hashtable_partialsums(hashfirstcodes.table,timer);
  hashfirstcodes.suftab = gt_malloc((size_t) psum *
                                    sizeof (*hashfirstcodes.suftab));
  hashfirstcodes.finalpsum = psum;
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "insert suffixes into suffix table",
                                        stdout);
  }
  getencseqkmers_twobitencoding(encseq,
                                GT_READMODE_FORWARD,
                                kmersize,
                                spmopt,
                                false,
                                gt_insertallcodeswithhashtable,
                                &hashfirstcodes,
                                NULL,
                                NULL);
  gt_free(hashfirstcodes.suftab);
  gt_uint64hashtable_delete(hashfirstcodes.table);
  if (timer != NULL)
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
}
