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

#ifndef ENCODEDSEQUENCE_H
#define ENCODEDSEQUENCE_H

#include "core/alphabet.h"
#include "core/chardef.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/symboldef.h"
#include "core/filelengthvalues.h"
#include "core/disc_distri.h"
#include "encodedsequence_api.h"
#include "core/seqpos.h"
#include "intcode-def.h"
#include "core/intbits.h"
#include "core/readmode.h"
#include "core/logger.h"

#ifdef SKDEBUG
#define GT_CHECKENCCHAR(CC,ENCSEQ,POS,READMODE)\
        {\
          GtUchar cctmp = gt_encodedsequence_getencodedchar(ENCSEQ,POS, \
                                                            READMODE);\
          if ((CC) != cctmp)\
          {\
            printf("file %s, line %d: pos = %lu:cc = %u != %u = ccreal\n",\
                   __FILE__,__LINE__,\
                   (unsigned long) (POS),\
                   (unsigned int) (CC),\
                   (unsigned int) cctmp);\
            exit(GT_EXIT_PROGRAMMING_ERROR);\
          }\
        }
#else
#define GT_CHECKENCCHAR(CC,ENCSEQ,POS,READMODE)
#endif

typedef struct
{
  Seqpos leftpos,
         rightpos;
} GtSequencerange;          /* \Typedef{Sequencerange} */

typedef struct
{
  Twobitencoding tbe;           /* two bit encoding */
  unsigned int unitsnotspecial; /* units which are not special */
  Seqpos position;
} EndofTwobitencoding;

typedef struct Specialrangeiterator Specialrangeiterator;

void extract2bitenc(bool fwd,
                    EndofTwobitencoding *ptbe,
                    const GtEncodedsequence *encseq,
                    GtEncodedsequenceScanstate *esr,
                    Seqpos startpos);

typedef struct
{
  unsigned int common;
  bool leftspecial, rightspecial;
  Seqpos finaldepth;
} GtCommonunits;

int compareTwobitencodings(bool fwd,
                           bool complement,
                           GtCommonunits *commonunits,
                           const EndofTwobitencoding *ptbe1,
                           const EndofTwobitencoding *ptbe2);

uint64_t detencseqofsatviatables(int kind,
                                 Seqpos totallength,
                                 unsigned long numofdbfiles,
                                 unsigned long lengthofdbfilenames,
                                 Seqpos specialranges,
                                 unsigned int numofchars);

void plainseq2bytecode(GtUchar *bytecode,const GtUchar *seq,unsigned long len);

void sequence2bytecode(GtUchar *dest,const GtEncodedsequence *encseq,
                       Seqpos startindex,Seqpos len);

int flushencseqfile(const GtStr *indexname,GtEncodedsequence *encseq,GtError*);

void checkallsequencedescriptions(const GtEncodedsequence *encseq);

GtEncodedsequence *plain2encodedsequence(bool withrange,
                                         const GtUchar *seq1,
                                         Seqpos len1,
                                         const GtUchar *seq2,
                                         unsigned long len2,
                                         const GtAlphabet *alpha,
                                         GtLogger *logger);

Specialrangeiterator *newspecialrangeiterator(const GtEncodedsequence *encseq,
                                              bool moveforward);

bool hasspecialranges(const GtEncodedsequence *encseq);

bool hasfastspecialrangeenumerator(const GtEncodedsequence *encseq);

bool possibletocmpbitwise(const GtEncodedsequence *encseq);

bool nextspecialrangeiterator(GtSequencerange *range,Specialrangeiterator *sri);

void freespecialrangeiterator(Specialrangeiterator **sri);

/*@null@*/ const char *encseqaccessname(const GtEncodedsequence *encseq);

Codetype extractprefixcode(unsigned int *unitsnotspecial,
                           const GtEncodedsequence *encseq,
                           const Codetype *filltable,
                           GtReadmode readmode,
                           GtEncodedsequenceScanstate *esr,
                           const Codetype **multimappower,
                           Seqpos frompos,
                           unsigned int prefixlength);

int comparewithonespecial(bool *leftspecial,
                          bool *rightspecial,
                          const GtEncodedsequence *encseq,
                          bool fwd,
                          bool complement,
                          Seqpos pos1,
                          Seqpos pos2,
                          Seqpos depth,
                          Seqpos maxdepth);

int compareEncseqsequences(GtCommonunits *commonunits,
                           const GtEncodedsequence *encseq,
                           bool fwd,
                           bool complement,
                           GtEncodedsequenceScanstate *esr1,
                           GtEncodedsequenceScanstate *esr2,
                           Seqpos pos1,
                           Seqpos pos2,
                           Seqpos depth);

int compareEncseqsequencesmaxdepth(GtCommonunits *commonunits,
                                   const GtEncodedsequence *encseq,
                                   bool fwd,
                                   bool complement,
                                   GtEncodedsequenceScanstate *esr1,
                                   GtEncodedsequenceScanstate *esr2,
                                   Seqpos pos1,
                                   Seqpos pos2,
                                   Seqpos depth,
                                   Seqpos maxdepth);

/* some check functions called in test-encseq.c */

int multicharactercompare(const GtEncodedsequence *encseq,
                          bool fwd,
                          bool complement,
                          GtEncodedsequenceScanstate *esr1,
                          Seqpos pos1,
                          GtEncodedsequenceScanstate *esr2,
                          Seqpos pos2);

void checkextractunitatpos(const GtEncodedsequence *encseq,
                           bool fwd,bool complement);

void checkextractspecialbits(const GtEncodedsequence *encseq,bool fwd);

void multicharactercompare_withtest(const GtEncodedsequence *encseq,
                                    bool fwd,
                                    bool complement,
                                    GtEncodedsequenceScanstate *esr1,
                                    Seqpos pos1,
                                    GtEncodedsequenceScanstate *esr2,
                                    Seqpos pos2);

void showsequenceatstartpos(FILE *fp,
                            bool fwd,
                            bool complement,
                            const GtEncodedsequence *encseq,
                            Seqpos startpos);

bool containsspecial(const GtEncodedsequence *encseq,
                     bool moveforward,
                     GtEncodedsequenceScanstate *esrspace,
                     Seqpos startpos,
                     Seqpos len);

int getsatforcevalue(const char *str,GtError *err);

/* check if the marked positions are correct */

void checkmarkpos(const GtEncodedsequence *encseq);

/* for a array of recordseparator, obtain the sequence
 * number from the given position */

unsigned long getrecordnumSeqpos(const Seqpos *recordseps,
                                 unsigned long numofrecords,
                                 Seqpos totalwidth,
                                 Seqpos position);

/* for a given GtEncodedsequence mapped with withssptab=true, obtain the
 * sequence number from the given position */

unsigned long getencseqfrompos2seqnum(const GtEncodedsequence *encseq,
                                      Seqpos position);

/* here are some functions to extract the different components of the
 * specialcharinfo included in encseq */

Seqpos getencseqspecialcharacters(const GtEncodedsequence *encseq);

Seqpos getencseqspecialranges(const GtEncodedsequence *encseq);

Seqpos getencseqrealspecialranges(const GtEncodedsequence *encseq);

Seqpos getencseqlengthofspecialprefix(const GtEncodedsequence *encseq);

Seqpos getencseqlengthofspecialsuffix(const GtEncodedsequence *encseq);

/* In case an GtEncodedsequence is not mapped, we still need to obtain the
   Specialcharainfo. This is done by the following function */

int readSpecialcharinfo(Specialcharinfo *specialcharinfo,
                        const GtStr *indexname,GtError *err);

/* some functions to obtain some components from the Alphabet pointed to
   by encseq->alpha */

unsigned int gt_encodedsequence_alphabetnumofchars(
                                               const GtEncodedsequence *encseq);

const GtUchar *gt_encodedsequence_alphabetsymbolmap(
                                               const GtEncodedsequence *encseq);

const GtUchar *gt_encodedsequence_alphabetcharacters(
                                               const GtEncodedsequence *encseq);

GtUchar gt_encodedsequence_alphabetwildcardshow(
                                               const GtEncodedsequence *encseq);

unsigned long getencseqcharactercount(const GtEncodedsequence *encseq,
                                      GtUchar cc);

/* some functions to remove reference from an GtEncodedsequence to prevent that
   the referenced alphabet or filenametab are freed */

void removealpharef(GtEncodedsequence *encseq);

void removefilenametabref(GtEncodedsequence *encseq);

void showgetencodedcharcounters(void);

void gt_showsequencefeatures(GtLogger *logger,
                             const GtEncodedsequence *encseq,
                             bool withfilenames);

unsigned long determinelengthofdbfilenames(const GtStrArray *filenametab);

FILE *opendestabfile(const GtStr *indexname,const char *mode,GtError *err);

FILE *openssptabfile(const GtStr *indexname,const char *mode,GtError *err);

int comparetwosuffixes(const GtEncodedsequence *encseq,
                       GtReadmode readmode,
                       Seqpos *maxlcp,
                       bool specialsareequal,
                       bool specialsareequalatdepth0,
                       Seqpos maxdepth,
                       Seqpos start1,
                       Seqpos start2,
                       GtEncodedsequenceScanstate *esr1,
                       GtEncodedsequenceScanstate *esr2);

int comparetwostrings(const GtEncodedsequence *encseq,
                      bool fwd,
                      bool complement,
                      Seqpos *maxcommon,
                      Seqpos pos1,
                      Seqpos pos2,
                      Seqpos maxdepth);

int comparetwostringsgeneric(const GtEncodedsequence *encseq,
                             bool fwd,
                             bool complement,
                             Seqpos *maxcommon,
                             Seqpos pos1,
                             Seqpos pos2,
                             Seqpos depth,
                             Seqpos maxdepth);

#endif
