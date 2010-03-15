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
#include "core/symboldef.h"
#include "core/logger.h"
#include "seqpos-def.h"
#include "intcode-def.h"

typedef struct
{
  Seqpos left;
  unsigned long nonspecialsinbucket,
                specialsinbucket,
                ordercode;
} Bucketspecification;

typedef struct Bcktab Bcktab;

Bcktab *mapbcktab(const GtStr *indexname,
                  unsigned int numofchars,
                  unsigned int prefixlength,
                  GtError *err);

void bcktab_delete(Bcktab **bcktab);

Bcktab *allocBcktab(unsigned int numofchars,
                    unsigned int prefixlength,
                    bool storespecialcodes,
                    GtLogger *logger,
                    GtError *err);

void updatebckspecials(Bcktab *bcktab,
                       Codetype code,
                       unsigned int numofchars,
                       unsigned int prefixindex);

Codetype codedownscale(const Bcktab *bcktab,
                       Codetype code,
                       unsigned int prefixindex,
                       unsigned int maxprefixlen);

void addfinalbckspecials(Bcktab *bcktab,unsigned int numofchars,
                         Seqpos specialcharacters);

int bcktab2file(FILE *fp,const Bcktab *bcktab,GtError *err);

unsigned int calcbucketboundsparts(Bucketspecification *bucketspec,
                                   const Bcktab *bcktab,
                                   Codetype code,
                                   Codetype maxcode,
                                   Seqpos totalwidth,
                                   unsigned int rightchar,
                                   unsigned int numofchars);

Seqpos calcbucketrightbounds(const Bcktab *bcktab,
                             Codetype code,
                             Codetype maxcode,
                             Seqpos totalwidth);

unsigned long distpfxidxpartialsums(const Bcktab *bcktab,Codetype code,
                                    unsigned int lowerbound);

void calcbucketboundaries(Bucketspecification *bucketspec,
                          const Bcktab *bcktab,
                          Codetype code);

void determinemaxbucketsize(Bcktab *bcktab,
                            const Codetype mincode,
                            const Codetype maxcode,
                            Seqpos partwidth,
                            unsigned int numofchars,
                            bool hashexceptions,
                            Seqpos totallength);/* relevant for hashexception */

void bcktab_showlog2info(const Bcktab *bcktab, GtLogger *logger);

unsigned int singletonmaxprefixindex(const Bcktab *bcktab,Codetype code);

unsigned long bcktab_specialsmaxbucketsize(const Bcktab *bcktab);

unsigned long bcktab_nonspecialsmaxbucketsize(const Bcktab *bcktab);

unsigned int bcktab_optimalnumofbits(unsigned short *logofremaining,
                                     const Bcktab *bcktab);

unsigned int pfxidx2lcpvalues(unsigned int *minprefixindex,
                              uint8_t *lcpsubtab,
                              unsigned long specialsinbucket,
                              const Bcktab *bcktab,
                              Codetype code);

const Codetype **bcktab_multimappower(const Bcktab *bcktab);

Codetype bcktab_filltable(const Bcktab *bcktab,unsigned int idx);

Seqpos *bcktab_leftborder(Bcktab *bcktab);

Codetype bcktab_numofallcodes(const Bcktab *bcktab);

uint64_t sizeofbuckettable(unsigned int numofchars,
                           unsigned int prefixlength);

unsigned int bcktab_prefixlength(const Bcktab *bcktab);

void bcktab_leftborderpartialsums(Bcktab *bcktab,Seqpos numofsuffixestosort);

#ifdef SKDEBUG
void checkcountspecialcodes(const Bcktab *bcktab);

void consistencyofsuffix(int line,
                         const Encodedsequence *encseq,
                         Readmode readmode,
                         const Bcktab *bcktab,
                         unsigned int numofchars,
                         const Suffixwithcode *suffix);
#endif
#endif
