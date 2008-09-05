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
#include "absdfstrans-def.h"

typedef struct Limdfsresources Limdfsresources;

Limdfsresources *newLimdfsresources(const void *genericindex,
                                    const Matchbound **mbtab,
                                    unsigned int maxdepth,
                                    const Encodedsequence *encseq,
                                    bool withesa,
                                    bool nowildcards,
                                    unsigned long maxintervalwidth,
                                    unsigned int mapsize,
                                    Seqpos totallength,
                                    unsigned long maxpathlength,
                                    Processmatch processmatch,
                                    void *processmatchinfo,
                                    Processresult processresult,
                                    void *patterninfo,
                                    const AbstractDfstransformer *adfst);

void freeLimdfsresources(Limdfsresources **ptrlimdfsresources,
                         const AbstractDfstransformer *adfst);

void indexbasedapproxpatternmatching(Limdfsresources *limdfsresources,
                                     bool rcmatch,
                                     const Uchar *pattern,
                                     unsigned long patternlength,
                                     unsigned long maxdistance,
                                     unsigned long maxintervalwidth,
                                     bool skpp,
                                     const AbstractDfstransformer *adfst);

unsigned long genericmstats(const Limdfsresources *limdfsresources,
                            const Uchar *qstart,
                            const Uchar *qend);

void indexbasedexactpatternmatching(const Limdfsresources *limdfsresources,
                                    bool rcmatch,
                                    const Uchar *pattern,
                                    unsigned long patternlength,
                                    Processmatch processmatch,
                                    void *processmatchinfo);

Seqpos bound2startpos(const Limdfsresources *limdfsresources,
                      Seqpos bound,unsigned long matchlength);

Uchar limdfsgetencodedchar(const Limdfsresources *limdfsresources,
                           Seqpos pos,
                           Readmode readmode);

Seqpos getlastbound(const Limdfsresources *limdfsresources,Seqpos rightbound);

bool intervalwidthleq(const Limdfsresources *limdfsresources,
                      Seqpos leftbound,Seqpos rightbound);

DECLAREARRAYSTRUCT(Seqpos);

ArraySeqpos *fromitv2sortedmatchpositions(Limdfsresources *limdfsresources,
                                          bool rcmatch,
                                          Seqpos leftbound,
                                          Seqpos rightbound,
                                          unsigned long offset);

#endif
