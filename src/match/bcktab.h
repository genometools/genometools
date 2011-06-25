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

#ifndef BCKTAB_H
#define BCKTAB_H

#include "core/error.h"
#include "core/intbits.h"
#include "core/codetype.h"

typedef struct
{
  unsigned long left,
                nonspecialsinbucket,
                specialsinbucket;
} GtBucketspecification;

typedef struct GtLeftborder GtLeftborder;

typedef struct GtBcktab GtBcktab;

void gt_bcktab_leftborder_addcode(GtLeftborder *lb,GtCodetype code);

unsigned long gt_bcktab_leftborder_insertionindex(GtLeftborder *lb,
                                                  GtCodetype code);

void gt_bcktab_leftborder_assign(GtLeftborder *lb,GtCodetype code,
                                 unsigned long value);

unsigned long gt_bcktab_get(const GtBcktab *bcktab,GtCodetype code);

GtBcktab *gt_bcktab_map(const char *indexname,
                        unsigned int numofchars,
                        unsigned int prefixlength,
                        bool withspecialsuffixes,
                        GtError *err);

void gt_bcktab_delete(GtBcktab *bcktab);

GtBcktab *gt_bcktab_alloc(unsigned int numofchars,
                          unsigned int prefixlength,
                          bool storespecialcodes,
                          bool withspecialsuffixes,
                          GtError *err);

void gt_bcktab_updatespecials(GtBcktab *bcktab,
                              GtCodetype code,
                              unsigned int numofchars,
                              unsigned int prefixindex);

GtCodetype gt_bcktab_codedownscale(const GtBcktab *bcktab,
                                   GtCodetype code,
                                   unsigned int prefixindex,
                                   unsigned int maxprefixlen);

void gt_bcktab_addfinalspecials(GtBcktab *bcktab,
                                unsigned int numofchars,
                                unsigned long specialcharacters);

int gt_bcktab_flush_to_file(FILE *fp,const GtBcktab *bcktab,GtError *err);

unsigned int gt_bcktab_calcboundsparts(GtBucketspecification *bucketspec,
                                       const GtBcktab *bcktab,
                                       GtCodetype code,
                                       GtCodetype maxcode,
                                       unsigned long totalwidth,
                                       unsigned int rightchar,
                                       unsigned int numofchars);

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
                                unsigned long partwidth,
                                unsigned int numofchars);

unsigned int gt_bcktab_singletonmaxprefixindex(const GtBcktab *bcktab,
                                               GtCodetype code);

unsigned int gt_bcktab_pfxidx2lcpvalues_uint8(unsigned int *minprefixindex,
                                              uint8_t *smalllcpvalues,
                                              unsigned long specialsinbucket,
                                              const GtBcktab *bcktab,
                                              GtCodetype code);

unsigned int gt_bcktab_pfxidx2lcpvalues_ulong(unsigned int *minprefixindex,
                                              unsigned long *bucketoflcpvalues,
                                              unsigned long specialsinbucket,
                                              const GtBcktab *bcktab,
                                              GtCodetype code);

const GtCodetype **gt_bcktab_multimappower(const GtBcktab *bcktab);

GtCodetype gt_bcktab_filltable(const GtBcktab *bcktab,unsigned int idx);

GtLeftborder *gt_bcktab_leftborder(GtBcktab *bcktab);

GtCodetype gt_bcktab_numofallcodes(const GtBcktab *bcktab);

uint64_t gt_bcktab_sizeoftable(unsigned int numofchars,
                               unsigned int prefixlength,
                               bool withspecialsuffixes);

unsigned long gt_bcktab_sizeofworkspace(unsigned int prefixlength);

unsigned int gt_bcktab_prefixlength(const GtBcktab *bcktab);

unsigned long gt_bcktab_leftborderpartialsums(
                             unsigned long *saved_bucketswithoutwholeleaf,
                             unsigned long *numofsuffixestosort,
                             GtBcktab *bcktab,
                             const GtBitsequence *markwholeleafbuckets);

size_t gt_bcktab_sizeforlcpvalues(const GtBcktab *bcktab);

void gt_bcktab_showleftborder(const GtBcktab *bcktab);

GtCodetype gt_bcktab_findfirstlarger(const GtBcktab *bcktab,
                                     unsigned long suftaboffset);

#ifdef SKDEBUG
void gt_bcktab_checkcountspecialcodes(const GtBcktab *bcktab);

void gt_bcktab_consistencyofsuffix(int line,
                                   const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   const GtBcktab *bcktab,
                                   unsigned int numofchars,
                                   const Suffixwithcode *suffix);
#endif
#endif
