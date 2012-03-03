/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_LCPVALUES_H
#define SFX_LCPVALUES_H
#include "core/error_api.h"
#include "core/codetype.h"
#include "core/unused_api.h"
#include "core/encseq.h"
#include "extended/rmq.h"
#include "lcpoverflow.h"
#include "bcktab.h"
#include "sfx-suffixgetset.h"

typedef struct GtOutlcpinfo GtOutlcpinfo;

typedef struct
{
#ifndef NDEBUG
  GtBitsequence *isset;
#endif
  GtLcpvaluetype *bucketoflcpvalues;
  unsigned long numofentries,
                numoflargelcpvalues,
                lcptaboffset; /* This can be positive when the lcp-values
                                 for an entire range of suffixes (covering
                                 more than one bucket) must
                                 be sorted, as in the case of the difference
                                 cover. In other cases, only the lcp-values
                                 for a single bucket must be computed, which
                                 means that the value is 0. */
} GtLcpvalues;

typedef void (*GtFinalProcessBucket)(void *,
                                     const GtSuffixsortspace *,
                                     const GtLcpvalues *,
                                     unsigned long,
                                     unsigned long,
                                     unsigned long);

/*@unused@*/ static inline void gt_lcptab_update(GtLcpvalues *tableoflcpvalues,
                                                 unsigned long subbucketleft,
                                                 unsigned long idx,
                                                 unsigned long value)
{
  gt_assert (tableoflcpvalues != NULL &&
             tableoflcpvalues->bucketoflcpvalues != NULL &&
             tableoflcpvalues->lcptaboffset+subbucketleft+idx <
             tableoflcpvalues->numofentries);
#ifndef NDEBUG
  if (tableoflcpvalues->isset != NULL)
  {
    GT_SETIBIT(tableoflcpvalues->isset,
               tableoflcpvalues->lcptaboffset+subbucketleft+idx);
  }
#endif
  gt_assert(value <= GT_LCPVALUE_MAX);
  tableoflcpvalues->bucketoflcpvalues[tableoflcpvalues->lcptaboffset +
                                      subbucketleft + idx]
                                      = (GtLcpvaluetype) value;
  if (value >= (unsigned long) LCPOVERFLOW)
  {
    tableoflcpvalues->numoflargelcpvalues++; /* this may overcount as there may
                                                be some value at index <idx>
                                                which was already overflowing */
  }
}

/*@unused@*/ static inline unsigned long gt_lcptab_getvalue(
                                        const GtLcpvalues *tableoflcpvalues,
                                        unsigned long subbucketleft,
                                        unsigned long idx)
{
  gt_assert (tableoflcpvalues != NULL &&
             tableoflcpvalues->bucketoflcpvalues != NULL &&
             tableoflcpvalues->lcptaboffset+subbucketleft+idx <
             tableoflcpvalues->numofentries);
  gt_assert(tableoflcpvalues->isset == NULL ||
            GT_ISIBITSET(tableoflcpvalues->isset,
                         tableoflcpvalues->lcptaboffset+subbucketleft+idx));
  return (unsigned long) tableoflcpvalues->bucketoflcpvalues
                           [tableoflcpvalues->lcptaboffset+subbucketleft+idx];
}

const GtLcpvaluetype *gt_lcptab_getptr(const GtLcpvalues *tableoflcpvalues,
                                       unsigned long subbucketleft);

GtOutlcpinfo *gt_Outlcpinfo_new(const char *indexname,
                                unsigned int numofchars,
                                unsigned int prefixlength,
                                bool withdistribution,
                                bool swallow_tail_lcpvalues,
                                GtFinalProcessBucket final_process_bucket,
                                void *final_process_bucket_info,
                                GtError *err);

size_t gt_Outlcpinfo_size(const GtOutlcpinfo *outlcpinfo);

void gt_Outlcpinfo_reinit(GtOutlcpinfo *outlcpinfo,
                          unsigned int numofchars,
                          unsigned int prefixlength,
                          unsigned long numoflcpvalues);

void gt_Outlcpinfo_delete(GtOutlcpinfo *outlcpinfo);

unsigned long gt_Outlcpinfo_numoflargelcpvalues(const GtOutlcpinfo *outlcpinfo);

double gt_Outlcpinfo_lcptabsum(const GtOutlcpinfo *outlcpinfo);

void gt_Outlcpinfo_numsuffixes2output_set(GtOutlcpinfo *outlcpinfo,
                                          unsigned long numsuffixes2output);

unsigned long gt_Outlcpinfo_maxbranchdepth(const GtOutlcpinfo *outlcpinfo);

void gt_Outlcpinfo_prebucket(GtOutlcpinfo *outlcpinfo,
                             GtCodetype code,
                             unsigned long lcptaboffset);

void gt_Outlcpinfo_nonspecialsbucket(GtOutlcpinfo *outlcpinfo,
                                     unsigned int prefixlength,
                                     const GtSuffixsortspace *sssp,
                                     GtLcpvalues *tableoflcpvalues,
                                     const GtBucketspecification *bucketspec,
                                     GtCodetype code);

void gt_Outlcpinfo_postbucket(GtOutlcpinfo *outlcpinfo,
                              unsigned int prefixlength,
                              const GtSuffixsortspace *sssp,
                              const GtBcktab *bcktab,
                              const GtBucketspecification *bucketspec,
                              GtCodetype code);

GtLcpvalues *gt_Outlcpinfo_resizereservoir(GtOutlcpinfo *outlcpinfo,
                                           const GtBcktab *bcktab);

GtLcpvalues *gt_Outlcpinfo_lcpvalues_ref(GtOutlcpinfo *outlcpinfo);

void gt_Outlcpinfo_check_lcpvalues(const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   const GtSuffixsortspace *sortedsample,
                                   unsigned long effectivesamplesize,
                                   const GtOutlcpinfo *outlcpinfosample,
                                   bool checkequality);

GtRMQ *gt_lcpvalues_rmq_new(const GtLcpvalues *samplelcpvalues);

#endif
