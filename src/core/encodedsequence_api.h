/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#ifndef ENCODEDSEQUENCE_API_H
#define ENCODEDSEQUENCE_API_H

#include "core/alphabet.h"
#include "core/chardef.h"
#include "core/logger.h"
#include "core/progress_timer.h"
#include "core/readmode.h"
#include "core/seqpos.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/symboldef.h"

typedef struct GtEncodedsequence GtEncodedsequence;
typedef struct GtEncodedsequenceScanstate GtEncodedsequenceScanstate;

#undef GT_INLINEDENCSEQ
#ifdef GT_INLINEDENCSEQ
#include "core/encodedsequence_rep.h"
#endif

/*@null@*/
GtEncodedsequence* gt_encodedsequence_new_from_files(
                                                  GtProgressTimer *sfxprogress,
                                                  const GtStr *str_indexname,
                                                  const GtStr *str_smap,
                                                  const GtStr *str_sat,
                                                  const GtStrArray *filenametab,
                                                  bool isdna,
                                                  bool isprotein,
                                                  bool isplain,
                                                  bool outtistab,
                                                  bool outdestab,
                                                  bool outsdstab,
                                                  bool outssptab,
                                                  GtLogger *logger,
                                                  GtError *err);

/*@null@*/
GtEncodedsequence* gt_encodedsequence_new_from_index(bool withrange,
                                                     const GtStr *indexname,
                                                     bool withtistab,
                                                     bool withdestab,
                                                     bool withsdstab,
                                                     bool withssptab,
                                                     GtLogger *logger,
                                                     GtError *err);

#ifdef GT_INLINEDENCSEQ
#define            gt_encodedsequence_total_length(ENCSEQ) \
                     ((ENCSEQ)->totallength)
#else
Seqpos             gt_encodedsequence_total_length(
                                               const GtEncodedsequence *encseq);
#endif

#ifdef GT_INLINEDENCSEQ
#define            gt_encodedsequence_num_of_sequences(ENCSEQ) \
                     ((ENCSEQ)->numofdbsequences)
#else
unsigned long      gt_encodedsequence_num_of_sequences(
                                               const GtEncodedsequence *encseq);
#endif

#define GT_REVERSEPOS(TOTALLENGTH,POS) \
          ((TOTALLENGTH) - 1 - (POS))

#ifdef GT_INLINEDENCSEQ
#define GT_MAKECOMPL(CC) \
          (ISSPECIAL(CC) ? (CC) : (GtUchar) 3 - (CC))
/*@unused@*/ static inline
GtUchar            gt_encodedsequence_getencodedchar(
                                                const GtEncodedsequence *encseq,
                                                Seqpos pos,
                                                GtReadmode readmode)
{
  return (readmode == GT_READMODE_FORWARD)
          ? encseq->plainseq[pos]
          : ((readmode == GT_READMODE_REVERSE)
            ? encseq->plainseq[GT_REVERSEPOS(encseq->totallength,pos)]
            : ((readmode == GT_READMODE_COMPL)
              ? GT_MAKECOMPL(encseq->plainseq[pos])
              : GT_MAKECOMPL(encseq->plainseq[
                           GT_REVERSEPOS(encseq->totallength,pos)])
              )
            )
         ;
}
#define            gt_encodedsequence_extractencodedchar(ENCSEQ,POS,RM) \
                     gt_encodedsequence_getencodedchar(ENCSEQ,POS,RM)
#else
GtUchar            gt_encodedsequence_getencodedchar(
                                                const GtEncodedsequence *encseq,
                                                Seqpos pos,
                                                GtReadmode readmode);
GtUchar            gt_encodedsequence_extractencodedchar(
                                                const GtEncodedsequence *encseq,
                                                Seqpos pos,
                                                GtReadmode readmode);
#endif

#ifdef GT_INLINEDENCSEQ
#define            gt_encodedsequence_getencodedcharnospecial(ENCSEQ,POS,RM) \
                     gt_encodedsequence_getencodedchar(ENCSEQ,POS,RM)
#else
GtUchar            gt_encodedsequence_getencodedcharnospecial(
                                                const GtEncodedsequence *encseq,
                                                Seqpos pos,
                                                GtReadmode readmode);
#endif

#ifdef GT_INLINEDENCSEQ
#define            gt_encodedsequence_sequentialgetencodedchar(ENCSEQ, \
                                                     ENCSEQSTATE,POS,READMODE) \
                     gt_encodedsequence_getencodedchar(ENCSEQ,POS,READMODE)
#else
GtUchar            gt_encodedsequence_sequentialgetencodedchar(
                                                const GtEncodedsequence *encseq,
                                                GtEncodedsequenceScanstate *esr,
                                                Seqpos pos,
                                                GtReadmode readmode);
#endif

void               gt_encodedsequence_extract_substring(
                                                const GtEncodedsequence *encseq,
                                                GtUchar *buffer,
                                                Seqpos frompos,
                                                Seqpos topos);

typedef struct
{
  Seqpos seqstartpos,  /* the position of the first character in the encseq */
         seqlength;    /* the length of the sequence */
} GtSeqinfo;

/* Fills the <seqinfo> struct for the <seqnum>-th sequence in the <encseq>. */
void               gt_encodedsequence_seqinfo(const GtEncodedsequence *encseq,
                                              GtSeqinfo *seqinfo,
                                              unsigned long seqnum);

/* Returns a pointer to the description of the <seqnum>-th sequence in the
   <encseq>. The length of the returned string is written to the
   location pointed at by <desclen>. */
const char*        gt_encodedsequence_description(
                                                const GtEncodedsequence *encseq,
                                                unsigned long *desclen,
                                                unsigned long seqnum);

/* Returns the <GtAlphabet> associated with <encseq>. */
const GtAlphabet*  gt_encodedsequence_alphabet(const GtEncodedsequence *encseq);

/* Returns a <GtStrArray> of the names of the original sequence files
   contained in <encseq>. */
const GtStrArray*  gt_encodedsequence_filenames(
                                               const GtEncodedsequence *encseq);

void               gt_encodedsequence_delete(GtEncodedsequence *encseq);

/* TODO document: needed for efficient sequential reading */
GtEncodedsequenceScanstate* gt_encodedsequence_scanstate_new(void);
void                        gt_encodedsequence_scanstate_init(
                                                GtEncodedsequenceScanstate *esr,
                                                const GtEncodedsequence *encseq,
                                                GtReadmode readmode,
                                                Seqpos startpos);
void                        gt_encodedsequence_scanstate_initgeneric(
                                                GtEncodedsequenceScanstate *esr,
                                                const GtEncodedsequence *encseq,
                                                bool moveforward,
                                                Seqpos startpos);
void                        gt_encodedsequence_scanstate_delete(
                                               GtEncodedsequenceScanstate *esr);

#endif
