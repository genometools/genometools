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

/* TODO: what does this struct contain and how can it be used? */
typedef struct
{
  GtTwobitencoding tbe;           /* two bit encoding */
  unsigned int unitsnotspecial;   /* units which are not special */
  unsigned long position;
} GtEndofTwobitencoding;

/* TODO: what does this struct contain and how can it be used? */
typedef struct
{
  unsigned int common;
  bool leftspecial, rightspecial;
  unsigned long finaldepth;
} GtCommonunits;

/* The <GtSpecialrangeiterator> type. */
typedef struct GtSpecialrangeiterator GtSpecialrangeiterator;

/* Create a new <GtSpecialrangeiterator> for <encseq>. */
GtSpecialrangeiterator *
     gt_specialrangeiterator_new(const GtEncodedsequence *encseq,
                                 bool moveforward);

/* Make <sri> supply the next special range <range>. Returns true if another
   range was returned, false otherwise. */
bool gt_specialrangeiterator_next(GtSpecialrangeiterator *sri,
                                  GtRange *range);

/* Delete <sri> and free associated memory. */
void gt_specialrangeiterator_delete(GtSpecialrangeiterator *sri);

/* TODO: please document me */
void gt_encodedsequence_extract2bitenc(bool fwd,
                                       GtEndofTwobitencoding *ptbe,
                                       const GtEncodedsequence *encseq,
                                       GtEncodedsequenceScanstate *esr,
                                       unsigned long startpos);

/* TODO: please document me */
int gt_encodedsequence_compare_twobitencodings(bool fwd,
                                            bool complement,
                                            GtCommonunits *commonunits,
                                            const GtEndofTwobitencoding *ptbe1,
                                            const GtEndofTwobitencoding *ptbe2);

/* TODO: please document me */
void gt_encodedsequence_plainseq2bytecode(GtUchar *bytecode,
                                          const GtUchar *seq,
                                          unsigned long len);

/* TODO: please document me */
void gt_encodedsequence_sequence2bytecode(GtUchar *dest,
                                          const GtEncodedsequence *encseq,
                                          unsigned long startindex,
                                          unsigned long len);

/* Similar to <gt_error_check()>, this function exits with an error message
   if <encseq> returns inconsistent descriptions when compared to the
   destab. */
void gt_encodedsequence_check_descriptions(const GtEncodedsequence *encseq);

/* Similar to <gt_error_check()>, this function exits with an error message
   if <encseq> returns inconsistent marked positions. */
void gt_encodedsequence_check_markpos(const GtEncodedsequence *encseq);

/* TODO: please document me */
int gt_encodedsequence_check_specialranges(const GtEncodedsequence *encseq);

/* TODO: please document me */
int gt_encodedsequence_check_consistency(const GtEncodedsequence *encseq,
                                         const GtStrArray *filenametab,
                                         GtReadmode readmode,
                                         unsigned long scantrials,
                                         unsigned long multicharcmptrials,
                                         GtError *err);

/* Returns true is <encseq> has special ranges, false otherwise. */
bool gt_encodedsequence_has_specialranges(const GtEncodedsequence *encseq);

/* TODO: please document me */
bool gt_encodedsequence_has_fast_specialrangeenumerator(
                                               const GtEncodedsequence *encseq);

/* TODO: please document me */
bool gt_encodedsequence_bitwise_cmp_ok(const GtEncodedsequence *encseq);

/* TODO: please document me */
/*@null@*/
const char* gt_encodedsequence_accessname(const GtEncodedsequence *encseq);

GtCodetype gt_encodedsequence_extractprefixcode(unsigned int *unitsnotspecial,
                                               const GtEncodedsequence *encseq,
                                               const GtCodetype *filltable,
                                               GtReadmode readmode,
                                               GtEncodedsequenceScanstate *esr,
                                               const GtCodetype **multimappower,
                                               unsigned long frompos,
                                               unsigned int prefixlength);

int        gt_encodedsequence_compare(const GtEncodedsequence *encseq,
                                      GtCommonunits *commonunits,
                                      bool fwd,
                                      bool complement,
                                      GtEncodedsequenceScanstate *esr1,
                                      GtEncodedsequenceScanstate *esr2,
                                      unsigned long pos1,
                                      unsigned long pos2,
                                      unsigned long depth);

int        gt_encodedsequence_compare_maxdepth(const GtEncodedsequence *encseq,
                                               GtCommonunits *commonunits,
                                               bool fwd,
                                               bool complement,
                                               GtEncodedsequenceScanstate *esr1,
                                               GtEncodedsequenceScanstate *esr2,
                                               unsigned long pos1,
                                               unsigned long pos2,
                                               unsigned long depth,
                                               unsigned long maxdepth);

bool       gt_encodedsequence_contains_special(const GtEncodedsequence *encseq,
                                           bool moveforward,
                                           GtEncodedsequenceScanstate *esrspace,
                                           unsigned long startpos,
                                           unsigned long len);

/* Returns the sequence number from the given <position> for an array of of
   SEPARATOR positions <recordseps>.  */
unsigned long gt_encodedsequence_sep2seqnum(const unsigned long *recordseps,
                                            unsigned long numofrecords,
                                            unsigned long totalwidth,
                                            unsigned long position);

/* Returns the sequence number from the given <position> for a given
   GtEncodedsequence <encseq> mapped with withssptab=true. */
unsigned long gt_encodedsequence_pos2seqnum(const GtEncodedsequence *encseq,
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

/* Obtains the number of characters in the Alphabet associated with <encseq>.
   This saves one function call for extracting the alphabet pointer
   from <encseq> (performance reasons). */
unsigned int gt_encodedsequence_alphabetnumofchars(
                                               const GtEncodedsequence *encseq);

/* Obtains the symbolmap from the Alphabet associated with <encseq>.
   This saves one function call for extracting the alphabet pointer
   from <encseq> (performance reasons). */
const GtUchar *gt_encodedsequence_alphabetsymbolmap(
                                               const GtEncodedsequence *encseq);

/* Obtains an ordered array of characters from the Alphabet associated with
   <encseq>. This saves one function call for extracting the alphabet pointer
   from <encseq> (performance reasons). */
const GtUchar *gt_encodedsequence_alphabetcharacters(
                                               const GtEncodedsequence *encseq);

/* Obtains the character used for displaying wildcards from the Alphabet
   associated with <encseq>. This saves one function call for extracting the
   alphabet pointer from <encseq> (performance reasons). */
GtUchar gt_encodedsequence_alphabetwildcardshow(
                                               const GtEncodedsequence *encseq);

/* Returns the number of times that <cc> occurs in the sequences in <encseq>. */
unsigned long gt_encodedsequence_charcount(const GtEncodedsequence *encseq,
                                           GtUchar cc);

void gt_encodedsequence_show_features(const GtEncodedsequence *encseq,
                                      GtLogger *logger,
                                      bool withfilenames);

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
