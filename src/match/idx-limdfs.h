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

#include "core/readmode.h"
#include "procmatch.h"
#include "core/logger.h"
#include "absdfstrans-def.h"
#include "core/intbits.h"
#include "match/eis-voiditf.h"
#include "match/sarr-def.h"

typedef struct Genericindex Genericindex;

void genericindex_delete(Genericindex *genericindex);

const GtEncseq *genericindex_getencseq(const Genericindex
                                                *genericindex);

const Suffixarray *genericindex_getsuffixarray(const Genericindex
                                                *genericindex);

Genericindex *genericindex_new(const char *indexname,
                               bool withesa,
                               bool withencseq,
                               bool withdestab,
                               bool withssptab,
                               int userdefinedmaxdepth,
                               GtLogger *logger,
                               GtError *err);

typedef struct Limdfsresources Limdfsresources;

Limdfsresources *gt_newLimdfsresources(const Genericindex *genericindex,
                                    bool nowildcards,
                                    GtUword maxintervalwidth,
                                    GtUword maxpathlength,
                                    bool keepexpandedonstack,
                                    ProcessIdxMatch processmatch,
                                    void *processmatchinfo,
                                    Processresult processresult,
                                    void *patterninfo,
                                    const AbstractDfstransformer *adfst);

void gt_freeLimdfsresources(Limdfsresources **ptrlimdfsresources,
                         const AbstractDfstransformer *adfst);

bool gt_indexbasedapproxpatternmatching(Limdfsresources *limdfsresources,
                                     const GtUchar *pattern,
                                     GtUword patternlength,
                                     GtUword maxdistance,
                                     GtUword maxintervalwidth,
                                     bool skpp,
                                     const AbstractDfstransformer *adfst);

void gt_indexbasedmstats(Limdfsresources *limdfsresources,
                      const GtUchar *pattern,
                      GtUword patternlength,
                      const AbstractDfstransformer *adfst);

void gt_indexbasedspacedseeds(Limdfsresources *limdfsresources,
                           const GtUchar *pattern,
                           GtBitsequence seedbitvector,
                           GtUword seedweight,
                           const AbstractDfstransformer *adfst);

void gt_indexbasedlocali(Limdfsresources *limdfsresources,
                      GtWord matchscore,
                      GtWord mismatchscore,
                      GtWord gapstart,
                      GtWord gapextend,
                      GtUword threshold,
                      const GtUchar *query,
                      GtUword querylength,
                      const AbstractDfstransformer *adfst);

GtUword genericmstats(const Limdfsresources *limdfsresources,
                            const GtUchar *qstart,
                            const GtUchar *qend);

bool gt_indexbasedexactpatternmatching(const Limdfsresources *limdfsresources,
                                    const GtUchar *pattern,
                                    GtUword patternlength);

GtUchar gt_limdfs_getencodedchar(const Limdfsresources *limdfsresources,
                              GtUword pos,
                              GtReadmode readmode);

bool gt_intervalwidthleq(const Limdfsresources *limdfsresources,
                      GtUword leftbound,GtUword rightbound);

GtArrayGtUword *gt_fromitv2sortedmatchpositions(
                                          Limdfsresources *limdfsresources,
                                          GtUword leftbound,
                                          GtUword rightbound,
                                          GtUword offset);

const FMindex *genericindex_get_packedindex(const Genericindex *genericindex);

GtUword genericindex_get_totallength(const Genericindex *genericindex);

GtUword gt_indexbased_exact_pattern_count(
                                              const Genericindex *genericindex,
                                              const GtUchar *pattern,
                                              GtUword patternlength);

#endif
