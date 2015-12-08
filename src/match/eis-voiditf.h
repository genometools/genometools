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

#ifndef EIS_VOIDITF_H
#define EIS_VOIDITF_H

#include "core/unused_api.h"
#include "core/encseq_api.h"

#include "splititv.h"
#include "procmatch.h"

/* The following type is just used to stronlgy type the functions
   using the FMindex (rather than using void pointers) */

typedef struct FMindex FMindex;

GtUword gt_bwtseqfirstmatch(const FMindex *fmindex, GtUword bound);

/** Iterator for positions defined by a lower and upper bound
 */

typedef struct Bwtseqpositioniterator Bwtseqpositioniterator;

Bwtseqpositioniterator *gt_Bwtseqpositioniterator_new(const FMindex *fmindex,
                                                      GtUword lowerbound,
                                                      GtUword upperbound);

bool gt_Bwtseqpositioniterator_next(GtUword *pos,
                                    Bwtseqpositioniterator *bspi);

bool gt_BwtseqpositionwithoutSEPiterator_next(GtUword *pos,
                                              Bwtseqpositioniterator *bspi);

GtUchar gt_bwtseqgetsymbol(GtUword bound,const FMindex *fmindex);

void gt_Bwtseqpositioniterator_delete(Bwtseqpositioniterator *bspi);

typedef struct BwtSeqpositionextractor BwtSeqpositionextractor;

BwtSeqpositionextractor *gt_newBwtSeqpositionextractor(
                                                     const FMindex *voidbwtseq,
                                                     GtUword upperbound);

GtUword gt_BwtSeqpositionextractor_extract(BwtSeqpositionextractor *bspex,
                                                 GtUword lowerbound);

void gt_freeBwtSeqpositionextractor(BwtSeqpositionextractor *bspex);

typedef struct Bwtseqcontextiterator Bwtseqcontextiterator;

Bwtseqcontextiterator *gt_Bwtseqcontextiterator_new(const FMindex *fmindex,
                                                    GtUword bound);

GtUchar gt_Bwtseqcontextiterator_next(GtUword *bound,
                                      Bwtseqcontextiterator *bsci);

void gt_Bwtseqcontextiterator_delete(Bwtseqcontextiterator *bsci);

void gt_bwtrangesplitwithoutspecial(GtArrayBoundswithchar *bwci,
                                    GtUword *rangeOccs,
                                    const FMindex *fmindex,
                                    GtUword lbound,
                                    GtUword ubound);

FMindex *gt_loadvoidBWTSeqForSA(const char *indexname,
                                bool withpckbt,
                                GtError *err);

void gt_deletevoidBWTSeq(FMindex *packedindex);

/* the parameter is const void *, as this is required by the other
   indexed based methods */

GtUword gt_voidpackedindexuniqueforward(const void *fmindex,
                                              GT_UNUSED GtUword offset,
                                              GT_UNUSED GtUword left,
                                              GT_UNUSED GtUword right,
                                              GT_UNUSED GtUword
                                                        *witnessposition,
                                              const GtUchar *qstart,
                                              const GtUchar *qend);

/* the parameter is const void *, as this is required by the other
   indexed based methods */

GtUword gt_voidpackedindexmstatsforward(const void *fmindex,
                                              GT_UNUSED GtUword offset,
                                              GT_UNUSED GtUword left,
                                              GT_UNUSED GtUword right,
                                              GtUword *witnessposition,
                                              const GtUchar *qstart,
                                              const GtUchar *qend);

bool gt_pck_exactpatternmatching(const FMindex *fmindex,
                                 const GtUchar *pattern,
                                 GtUword patternlength,
                                 GtUword totallength,
                                 const GtUchar *dbsubstring,
                                 ProcessIdxMatch processmatch,
                                 void *processmatchinfo);

GtUword gt_voidpackedfindfirstmatchconvert(const FMindex *fmindex,
                                                 GtUword witnessbound,
                                                 GtUword matchlength);

typedef struct
{
  GtUword lowerbound, upperbound;
} Mbtab;

GtUword gt_bwtrangesplitallwithoutspecial(Mbtab *mbtab,
                                                GtUword *rangeOccs,
                                                const FMindex *fmindex,
                                                GtUword lbound,
                                                GtUword ubound);

unsigned int gt_bwtseq2maxdepth(const FMindex *fmindex);

unsigned int gt_bwtseq2numofchars(const FMindex *fmindex);

const Mbtab **gt_bwtseq2mbtab(const FMindex *fmindex);

/* this does currently not work, only for the root interval. This is due to the
 * sorting of the special chars, might be changed in future.
 * Only reliable information:
 * - Child intervals of standard alphabet
 * - type of special and the number, but not the interval (child intervals in
 *   special ranges do not necessarily form continuous intervals!) */
GtUword gt_bwtrangesplitallwithspecial(Mbtab *mbtab,
                                             GtUword *rangeOccs,
                                             const FMindex *voidBwtSeq,
                                             GtUword lbound,
                                             GtUword ubound);

GtUword gt_voidpackedindex_totallength_get(const FMindex *fmindex);

/* returns the length of the matching prefix +1, that is it returns the
 * length of the shortest absent prefix */
GtUword gt_pck_getShuStringLength(const FMindex* bwtSubject,
                                       const GtUchar* query,
                                       GtUword queryLength);

/* returns the ½gc-content of the subject assuming fwd and rev sequence in the
 * index so that |g| == |c| */
double gt_pck_getGCcontent(const FMindex *bwtSubject,
                           const GtAlphabet *alphabet);

/* returns the count of all non special chars in the index, corresponding to the
 * first row in which a special char prefixes the rotation */
GtUword gt_pck_get_nonspecial_count(const FMindex *index);

/* returns occurrences of special chars in the intervals of non special chars,
 * that is the number of rows that would be extended with a special. */
GtUword gt_pck_special_occ_in_nonspecial_intervals(const FMindex *index);

/* counts the exact occurrences of pattern in index returns 0 if pattern is not
 found */
GtUword gt_pck_exact_pattern_count(const FMindex *index,
                                         const GtUchar *pattern,
                                         GtUword patternlength);
#endif
