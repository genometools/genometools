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

#include "core/assert_api.h"
#include "core/error.h"
#include "core/str.h"
#include "core/types_api.h"
#include "core/logger.h"
#include "core/codetype.h"

typedef struct
{
  unsigned long left,
                nonspecialsinbucket,
                specialsinbucket;
} Bucketspecification;

typedef struct Bcktab Bcktab;

Bcktab *gt_mapbcktab(const char *indexname,
                  unsigned int numofchars,
                  unsigned int prefixlength,
                  GtError *err);

void gt_bcktab_delete(Bcktab *bcktab);

Bcktab *gt_allocBcktab(unsigned int numofchars,
                    unsigned int prefixlength,
                    bool storespecialcodes,
                    GtError *err);

void gt_updatebckspecials(Bcktab *bcktab,
                       GtCodetype code,
                       unsigned int numofchars,
                       unsigned int prefixindex);

GtCodetype gt_codedownscale(const Bcktab *bcktab,
                       GtCodetype code,
                       unsigned int prefixindex,
                       unsigned int maxprefixlen);

void gt_addfinalbckspecials(Bcktab *bcktab,unsigned int numofchars,
                         unsigned long specialcharacters);

int gt_bcktab2file(FILE *fp,const Bcktab *bcktab,GtError *err);

unsigned int gt_calcbucketboundsparts(Bucketspecification *bucketspec,
                                   const Bcktab *bcktab,
                                   GtCodetype code,
                                   GtCodetype maxcode,
                                   unsigned long totalwidth,
                                   unsigned int rightchar,
                                   unsigned int numofchars);

unsigned long gt_calcbucketrightbounds(const Bcktab *bcktab,
                             GtCodetype code,
                             GtCodetype maxcode,
                             unsigned long totalwidth);

unsigned long gt_distpfxidxpartialsums(const Bcktab *bcktab,GtCodetype code,
                                    unsigned int lowerbound);

void gt_calcbucketboundaries(Bucketspecification *bucketspec,
                          const Bcktab *bcktab,
                          GtCodetype code);

void gt_determinemaxbucketsize(Bcktab *bcktab,
                            const GtCodetype mincode,
                            const GtCodetype maxcode,
                            unsigned long partwidth,
                            unsigned int numofchars);

unsigned int gt_singletonmaxprefixindex(const Bcktab *bcktab,GtCodetype code);

unsigned long gt_bcktab_nonspecialsmaxbucketsize(const Bcktab *bcktab);

unsigned int gt_pfxidx2lcpvalues_uint8(unsigned int *minprefixindex,
                                       uint8_t *smalllcpvalues,
                                       unsigned long specialsinbucket,
                                       const Bcktab *bcktab,
                                       GtCodetype code);

unsigned int gt_pfxidx2lcpvalues_ulong(unsigned int *minprefixindex,
                                       unsigned long *bucketoflcpvalues,
                                       unsigned long specialsinbucket,
                                       const Bcktab *bcktab,
                                       GtCodetype code);

const GtCodetype **gt_bcktab_multimappower(const Bcktab *bcktab);

GtCodetype gt_bcktab_filltable(const Bcktab *bcktab,unsigned int idx);

unsigned long *gt_bcktab_leftborder(Bcktab *bcktab);

GtCodetype gt_bcktab_numofallcodes(const Bcktab *bcktab);

uint64_t gt_sizeofbuckettable(unsigned int numofchars,
                              unsigned int prefixlength);

unsigned long gt_sizeofbucketworkspace(unsigned int prefixlength);

unsigned int gt_bcktab_prefixlength(const Bcktab *bcktab);

unsigned long gt_bcktab_emptybuckets(const Bcktab *bcktab);

unsigned long gt_bcktab_leftborderpartialsums(Bcktab *bcktab);

size_t gt_bcktab_sizeforlcpvalues(const Bcktab *bcktab);

#ifdef SKDEBUG
void checkcountspecialcodes(const Bcktab *bcktab);

void consistencyofsuffix(int line,
                         const GtEncseq *encseq,
                         GtReadmode readmode,
                         const Bcktab *bcktab,
                         unsigned int numofchars,
                         const Suffixwithcode *suffix);
#endif
#endif
