/*
  Copyright (c) 2007-2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef BCKTAB_H
#define BCKTAB_H

#include <limits.h>
#include "core/error.h"
#include "core/intbits.h"
#include "core/codetype.h"
#include "core/logger_api.h"
#include "sfx-maprange.h"

#ifdef GT_LONGLCPVALUES
typedef unsigned long GtLcpvaluetype;
#define GT_LCPVALUE_MAX ULONG_MAX
#else
typedef uint32_t GtLcpvaluetype;
#define GT_LCPVALUE_MAX UINT32_MAX
#endif

typedef struct
{
  unsigned long left,
                nonspecialsinbucket,
                specialsinbucket;
} GtBucketspecification;

typedef struct
{
  uint32_t *uintbounds;
  unsigned long *ulongbounds;
} GtLeftborder;

typedef struct GtBcktab GtBcktab;

/*@unused@*/ static inline void gt_bcktab_leftborder_addcode(GtLeftborder *lb,
                                                             GtCodetype code)
{
  gt_assert(lb != NULL);
  if (lb->ulongbounds != NULL)
  {
    lb->ulongbounds[code]++;
  } else
  {
    gt_assert(lb->uintbounds[code] < (uint32_t) UINT_MAX);
    lb->uintbounds[code]++;
  }
}

/*@unused@*/ static inline unsigned long gt_bcktab_leftborder_insertionindex(
                                                  GtLeftborder *lb,
                                                  GtCodetype code)
{
  gt_assert(lb != NULL);
  if (lb->ulongbounds != NULL)
  {
    return --lb->ulongbounds[code];
  }
  gt_assert(lb->uintbounds != NULL);
  return (unsigned long) --lb->uintbounds[code];
}

void gt_bcktab_leftborder_assign(GtLeftborder *lb,GtCodetype code,
                                 unsigned long value);

unsigned long gt_bcktab_get_leftborder(const GtBcktab *bcktab,GtCodetype code);

GtBcktab *gt_bcktab_map(const char *indexname,
                        unsigned int numofchars,
                        unsigned int prefixlength,
                        unsigned long maxvalue,
                        bool withspecialsuffixes,
                        GtError *err);

void gt_bcktab_assignboundsforpart(GtBcktab *bcktab,
                                   GtCodetype mincode,
                                   GtCodetype maxcode);

int gt_bcktab_remap_all(GtBcktab *bcktab,GtError *err);

void gt_bcktab_delete(GtBcktab *bcktab);

GtBcktab *gt_bcktab_new(unsigned int numofchars,
                        unsigned int prefixlength,
                        unsigned long maxvalue,
                        bool storespecialcodes,
                        bool withspecialsuffixes,
                        GtLogger *logger,
                        GtError *err);

void gt_bcktab_updatespecials(GtBcktab *bcktab,
                              GtCodetype code,
                              unsigned int prefixindex);

GtCodetype gt_bcktab_codedownscale(const GtBcktab *bcktab,
                                   GtCodetype code,
                                   unsigned int prefixindex,
                                   unsigned int maxprefixlen);

int gt_bcktab_flush_to_file(FILE *fp,const GtBcktab *bcktab,GtError *err);

void gt_bcktab_storetmp(GtBcktab *bcktab);

void gt_bcktab_maprange_lb_cs(GtSfxmappedrangelist *sfxmrlist,GtBcktab *bcktab);

unsigned int gt_bcktab_calcboundsparts(GtBucketspecification *bucketspec,
                                       const GtBcktab *bcktab,
                                       GtCodetype code,
                                       GtCodetype maxcode,
                                       unsigned long totalwidth,
                                       unsigned int rightchar);

unsigned long gt_bcktab_calcrightbounds(const GtBcktab *bcktab,
                                        GtCodetype code,
                                        GtCodetype maxcode,
                                        unsigned long totalwidth);

unsigned long gt_bcktab_distpfxidxpartialsums(const GtBcktab *bcktab,
                                              GtCodetype code,
                                              unsigned int lowerbound);

void gt_bcktab_calcboundaries(GtBucketspecification *bucketspec,
                              const GtBcktab *bcktab,
                              GtCodetype code);

void gt_bcktab_determinemaxsize(GtBcktab *bcktab,
                                const GtCodetype mincode,
                                const GtCodetype maxcode,
                                unsigned long partwidth);

unsigned int gt_bcktab_singletonmaxprefixindex(const GtBcktab *bcktab,
                                               GtCodetype code);

unsigned int gt_bcktab_pfxidx2lcpvalues_uint8(unsigned int *minprefixindex,
                                              uint8_t *smalllcpvalues,
                                              unsigned long specialsinbucket,
                                              const GtBcktab *bcktab,
                                              GtCodetype code);

unsigned int gt_bcktab_pfxidx2lcpvalues_Lcpvaluetype(
                                              unsigned int *minprefixindex,
                                              GtLcpvaluetype *bucketoflcpvalues,
                                              unsigned long specialsinbucket,
                                              const GtBcktab *bcktab,
                                              GtCodetype code);

const GtCodetype **gt_bcktab_multimappower(const GtBcktab *bcktab);

GtCodetype gt_bcktab_filltable(const GtBcktab *bcktab,unsigned int idx);

GtLeftborder *gt_bcktab_leftborder(GtBcktab *bcktab);

GtCodetype gt_bcktab_numofallcodes(const GtBcktab *bcktab);

uint64_t gt_bcktab_sizeoftable(unsigned int numofchars,
                               unsigned int prefixlength,
                               unsigned long maxvalue,
                               bool withspecialsuffixes);

unsigned long gt_bcktab_sizeofworkspace(unsigned int prefixlength);

unsigned int gt_bcktab_prefixlength(const GtBcktab *bcktab);

unsigned long gt_bcktab_leftborderpartialsums(
                             unsigned long *saved_bucketswithoutwholeleaf,
                             unsigned long *numofsuffixestosort,
                             GtBcktab *bcktab);

size_t gt_bcktab_sizeforlcpvalues(const GtBcktab *bcktab);

unsigned long gt_bcktab_maxbucketsize(const GtBcktab *bcktab);

unsigned long gt_bcktab_nonspecialsmaxsize(const GtBcktab *bcktab);

void gt_bcktab_leftborder_show(const GtBcktab *bcktab);

GtCodetype gt_bcktab_findfirstlarger(const GtBcktab *bcktab,
                                     unsigned long suftaboffset);

#ifdef SKDEBUG
void gt_bcktab_checkcountspecialcodes(const GtBcktab *bcktab);

void gt_bcktab_consistencyofsuffix(int line,
                                   const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   const GtBcktab *bcktab,
                                   const Suffixwithcode *suffix);
#endif
#endif
