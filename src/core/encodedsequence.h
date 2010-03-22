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
#include "core/encodedsequence_api.h"
#include "core/intbits.h"
#include "core/logger.h"
#include "core/range.h"
#include "core/readmode.h"
#include "match/intcode-def.h"

typedef struct
{
  GtTwobitencoding tbe;           /* two bit encoding */
  unsigned int unitsnotspecial;   /* units which are not special */
  unsigned long position;
} GtEndofTwobitencoding;

typedef struct
{
  unsigned int common;
  bool leftspecial, rightspecial;
  unsigned long finaldepth;
} GtCommonunits;

/* GtSpecialrangeiterator */

typedef struct GtSpecialrangeiterator GtSpecialrangeiterator;

GtSpecialrangeiterator *
     gt_specialrangeiterator_new(const GtEncodedsequence *encseq,
                                 bool moveforward);

bool gt_specialrangeiterator_next(GtSpecialrangeiterator *sri,
                                  GtRange *range);

void gt_specialrangeiterator_delete(GtSpecialrangeiterator *sri);

/* XXX: clean up interface funcs below */

void gt_encodedsequence_extract2bitenc(bool fwd,
                                       GtEndofTwobitencoding *ptbe,
                                       const GtEncodedsequence *encseq,
                                       GtEncodedsequenceScanstate *esr,
                                       unsigned long startpos);

int gt_encodedsequence_compare_twobitencodings(bool fwd,
                                            bool complement,
                                            GtCommonunits *commonunits,
                                            const GtEndofTwobitencoding *ptbe1,
                                            const GtEndofTwobitencoding *ptbe2);

void gt_encodedsequence_plainseq2bytecode(GtUchar *bytecode,
                                          const GtUchar *seq,
                                          unsigned long len);

void gt_encodedsequence_sequence2bytecode(GtUchar *dest,
                                          const GtEncodedsequence *encseq,
                                          unsigned long startindex,
                                          unsigned long len);

void gt_encodedsequence_check_descriptions(const GtEncodedsequence *encseq);

bool gt_encodedsequence_has_specialranges(const GtEncodedsequence *encseq);

bool gt_encodedsequence_has_fast_specialrangeenumerator(
                                               const GtEncodedsequence *encseq);

bool gt_encodedsequence_bitwise_cmp_ok(const GtEncodedsequence *encseq);

/*@null@*/
const char* gt_encodedsequence_accessname(const GtEncodedsequence *encseq);

Codetype extractprefixcode(unsigned int *unitsnotspecial,
                           const GtEncodedsequence *encseq,
                           const Codetype *filltable,
                           GtReadmode readmode,
                           GtEncodedsequenceScanstate *esr,
                           const Codetype **multimappower,
                           unsigned long frompos,
                           unsigned int prefixlength);

int comparewithonespecial(bool *leftspecial,
                          bool *rightspecial,
                          const GtEncodedsequence *encseq,
                          bool fwd,
                          bool complement,
                          unsigned long pos1,
                          unsigned long pos2,
                          unsigned long depth,
                          unsigned long maxdepth);

int compareEncseqsequences(GtCommonunits *commonunits,
                           const GtEncodedsequence *encseq,
                           bool fwd,
                           bool complement,
                           GtEncodedsequenceScanstate *esr1,
                           GtEncodedsequenceScanstate *esr2,
                           unsigned long pos1,
                           unsigned long pos2,
                           unsigned long depth);

int compareEncseqsequencesmaxdepth(GtCommonunits *commonunits,
                                   const GtEncodedsequence *encseq,
                                   bool fwd,
                                   bool complement,
                                   GtEncodedsequenceScanstate *esr1,
                                   GtEncodedsequenceScanstate *esr2,
                                   unsigned long pos1,
                                   unsigned long pos2,
                                   unsigned long depth,
                                   unsigned long maxdepth);

/* some check functions called in test-encseq.c */

int multicharactercompare(const GtEncodedsequence *encseq,
                          bool fwd,
                          bool complement,
                          GtEncodedsequenceScanstate *esr1,
                          unsigned long pos1,
                          GtEncodedsequenceScanstate *esr2,
                          unsigned long pos2);

void checkextractunitatpos(const GtEncodedsequence *encseq,
                           bool fwd,bool complement);

void checkextractspecialbits(const GtEncodedsequence *encseq,bool fwd);

void multicharactercompare_withtest(const GtEncodedsequence *encseq,
                                    bool fwd,
                                    bool complement,
                                    GtEncodedsequenceScanstate *esr1,
                                    unsigned long pos1,
                                    GtEncodedsequenceScanstate *esr2,
                                    unsigned long pos2);

void showsequenceatstartpos(FILE *fp,
                            bool fwd,
                            bool complement,
                            const GtEncodedsequence *encseq,
                            unsigned long startpos);

bool containsspecial(const GtEncodedsequence *encseq,
                     bool moveforward,
                     GtEncodedsequenceScanstate *esrspace,
                     unsigned long startpos,
                     unsigned long len);

int getsatforcevalue(const char *str,GtError *err);

/* check if the marked positions are correct */

void checkmarkpos(const GtEncodedsequence *encseq);

/* for a array of recordseparator, obtain the sequence
 * number from the given position */

unsigned long getrecordnumSeqpos(const unsigned long *recordseps,
                                 unsigned long numofrecords,
                                 unsigned long totalwidth,
                                 unsigned long position);

/* for a given GtEncodedsequence mapped with withssptab=true, obtain the
 * sequence number from the given position */

unsigned long getencseqfrompos2seqnum(const GtEncodedsequence *encseq,
                                      unsigned long position);

/* here are some functions to extract the different components of the
 * specialcharinfo included in encseq */

unsigned long getencseqspecialcharacters(const GtEncodedsequence *encseq);

unsigned long getencseqspecialranges(const GtEncodedsequence *encseq);

unsigned long getencseqrealspecialranges(const GtEncodedsequence *encseq);

unsigned long getencseqlengthofspecialprefix(const GtEncodedsequence *encseq);

unsigned long getencseqlengthofspecialsuffix(const GtEncodedsequence *encseq);

/* In case an GtEncodedsequence is not mapped, we still need to obtain the
   Specialcharainfo. This is done by the following function */

int readGtSpecialcharinfo(GtSpecialcharinfo *specialcharinfo,
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
                       unsigned long *maxlcp,
                       bool specialsareequal,
                       bool specialsareequalatdepth0,
                       unsigned long maxdepth,
                       unsigned long start1,
                       unsigned long start2,
                       GtEncodedsequenceScanstate *esr1,
                       GtEncodedsequenceScanstate *esr2);

int comparetwostrings(const GtEncodedsequence *encseq,
                      bool fwd,
                      bool complement,
                      unsigned long *maxcommon,
                      unsigned long pos1,
                      unsigned long pos2,
                      unsigned long maxdepth);

int comparetwostringsgeneric(const GtEncodedsequence *encseq,
                             bool fwd,
                             bool complement,
                             unsigned long *maxcommon,
                             unsigned long pos1,
                             unsigned long pos2,
                             unsigned long depth,
                             unsigned long maxdepth);

#endif
