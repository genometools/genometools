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

#ifndef IDX_LIMDFS_H
#define IDX_LIMDFS_H

#include "core/arraydef.h"
#include "seqpos-def.h"
#include "readmode-def.h"
#include "procmatch.h"
#include "verbose-def.h"
#include "absdfstrans-def.h"

typedef struct Genericindex Genericindex;

void genericindex_delete(Genericindex *genericindex);

const Encodedsequence *genericindex_getencseq(const Genericindex
                                              *genericindex);

Genericindex *genericindex_new(const GtStr *indexname,
                               bool withesa,
                               bool withencseq,
                               int userdefinedmaxdepth,
                               Verboseinfo *verboseinfo,
                               GtError *err);

typedef struct Limdfsresources Limdfsresources;

Limdfsresources *newLimdfsresources(const Genericindex *genericindex,
                                    bool nowildcards,
                                    unsigned long maxintervalwidth,
                                    unsigned long maxpathlength,
                                    Processmatch processmatch,
                                    void *processmatchinfo,
                                    Processresult processresult,
                                    void *patterninfo,
                                    const AbstractDfstransformer *adfst);

void freeLimdfsresources(Limdfsresources **ptrlimdfsresources,
                         const AbstractDfstransformer *adfst);

bool indexbasedapproxpatternmatching(Limdfsresources *limdfsresources,
                                     const Uchar *pattern,
                                     unsigned long patternlength,
                                     unsigned long maxdistance,
                                     unsigned long maxintervalwidth,
                                     bool skpp,
                                     const AbstractDfstransformer *adfst);

void indexbasedmstats(Limdfsresources *limdfsresources,
                      const Uchar *pattern,
                      unsigned long patternlength,
                      const AbstractDfstransformer *adfst);

void indexbasedspacedseeds(Limdfsresources *limdfsresources,
                           const Uchar *pattern,
                           Bitstring seedbitvector,
                           unsigned long seedweight,
                           const AbstractDfstransformer *adfst);

void indexbasedlocali(Limdfsresources *limdfsresources,
                      long matchscore,
                      long mismatchscore,
                      long gapstart,
                      long gapextend,
                      unsigned long threshold,
                      const Uchar *query,
                      unsigned long querylength,
                      const AbstractDfstransformer *adfst);

unsigned long genericmstats(const Limdfsresources *limdfsresources,
                            const Uchar *qstart,
                            const Uchar *qend);

bool indexbasedexactpatternmatching(const Limdfsresources *limdfsresources,
                                    const Uchar *pattern,
                                    unsigned long patternlength);

Seqpos bound2startpos(const Limdfsresources *limdfsresources,
                      Seqpos bound,unsigned long matchlength);

Uchar limdfsgetencodedchar(const Limdfsresources *limdfsresources,
                           Seqpos pos,
                           Readmode readmode);

bool intervalwidthleq(const Limdfsresources *limdfsresources,
                      Seqpos leftbound,Seqpos rightbound);

ArraySeqpos *fromitv2sortedmatchpositions(Limdfsresources *limdfsresources,
                                          Seqpos leftbound,
                                          Seqpos rightbound,
                                          unsigned long offset);

#endif
