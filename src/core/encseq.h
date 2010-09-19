/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2010 Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
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

#ifndef ENCSEQ_H
#define ENCSEQ_H

#include "core/alphabet.h"
#include "core/chardef.h"
#include "core/codetype.h"
#include "core/disc_distri.h"
#include "core/encseq_api.h"
#include "core/encseq_access_type.h"
#include "core/filelengthvalues.h"
#include "core/intbits.h"
#include "core/range.h"
#include "core/readmode.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/types_api.h"

#define GT_REVERSEPOS(TOTALLENGTH,POS) \
          ((TOTALLENGTH) - 1 - (POS))

/* The following type stores a two bit encoding in <tbe> with information
  about the number of two bit units which do not store a special
  character in <unitsnotspecial>. To allow the comparison of these
  structures, the <position> in the sequence at which the encoding was
  extracted is also stored */
typedef struct
{
  GtTwobitencoding tbe;           /* two bit encoding */
  unsigned int unitsnotspecial;   /* units which are not special */
  unsigned long position;
} GtEndofTwobitencoding;

/* The following type stores the result of comparing a pair of twobit
  encodings. <common> stores the number of units which are common
  (either from the beginning or from the end. <leftspecial> is true
  iff the first twobit encoding contains a special character in the units
  it covers. <rightspecial> is true iff the second twobit encoding contains
  a special character in the units it covers. <finaldepth> stores
  result of the comparison, i.e. the longest common prefix of the
  comparison. */
typedef struct
{
  unsigned int common;
  bool leftspecial,
       rightspecial;
  unsigned long finaldepth;
} GtCommonunits;

/* The <GtSpecialrangeiterator> type. */
typedef struct GtSpecialrangeiterator GtSpecialrangeiterator;

/* Create a new <GtSpecialrangeiterator> for <encseq>. */
GtSpecialrangeiterator* gt_specialrangeiterator_new(const GtEncseq *encseq,
                                                    bool moveforward);

/* Make <sri> supply the next special range <range>. Returns true if another
  range was returned, false otherwise. */
bool gt_specialrangeiterator_next(GtSpecialrangeiterator *sri,
                                  GtRange *range);

/* Delete <sri> and free associated memory. */
void gt_specialrangeiterator_delete(GtSpecialrangeiterator *sri);

/* Returns the encoded representation of the character at position <pos> of
  <encseq> read in the direction as indicated by <readmode>.
  The function only works for the case that encodesequence[pos] does not
  contain a special character. */
GtUchar gt_encseq_get_encoded_char_nospecial(const GtEncseq *encseq,
                                             unsigned long pos,
                                             GtReadmode readmode);

/* The following function extracts from the byte sequence of length <len>
  pointed to by <seq> the sequence of len/4 bytecodes each
  encoding four consecutive characters in one byte. */
void gt_encseq_plainseq2bytecode(GtUchar *bytecode,
                                 const GtUchar *seq,
                                 unsigned long len);

/* The following function extracts from an encoded sequence a substring of
  length <len> beginning at position <startindex> and stores the result
  in the byte sequence encoding four consecutive characters in one byte. */
void gt_encseq_sequence2bytecode(GtUchar *dest,
                                 const GtEncseq *encseq,
                                 unsigned long startindex,
                                 unsigned long len);

/* Returns true is <encseq> has special ranges, false otherwise. */
bool gt_encseq_has_specialranges(const GtEncseq *encseq);

/* Return true if the representation of the <encseq> is based on two
  bit encoding */
bool gt_encseq_bitwise_cmp_ok(const GtEncseq *encseq);

/* Return the integer code of the sequence of length <prefixlength>
  beginning at position <frompos> represetned by <encseq>. <esr> is
  used for efficiently scanning the sequence of symbols. <filltable>
  and <multimappower> are used for completing the integer code
  in cases where prefix of length <prefixlength> contains a special
  character. */
GtCodetype gt_encseq_extractprefixcode(unsigned int *unitsnotspecial,
                                       const GtEncseq *encseq,
                                       const GtCodetype *filltable,
                                       GtReadmode readmode,
                                       GtEncseqReader *esr,
                                       const GtCodetype **multimappower,
                                       unsigned long frompos,
                                       unsigned int prefixlength);

/* The following function compares two substrings beginning
  at position <pos1>+<depth> and <pos2>+<depth> in <encseq>. 
  <esr1> and <esr2> refer to memory areas for storeing a GtEncseqReader.
  The information about the length
  of the longest common prefix is stored in <commonunits>. <fwd> and
  <complement> specify if the sequence is scanned in forward direction
  and if the complement of the sequence is to be considered. The
  return value is -1, 0 or 1 depending on whether the sequence beginning at
  position <pos1>+<depth> is smaller than, equal to, or larger than the 
  sequence beginning at position <pos2>+<depth>. If <madepth> is 0, then the 
  entire suffixes are compared. Otherwise, the comparison is restricted to 
  the prefixes of length <maxdepth>. */
int gt_encseq_compare_viatwobitencoding(const GtEncseq *encseq,
                                        GtCommonunits *commonunits,
                                        GtReadmode readmode,
                                        GtEncseqReader *esr1,
                                        GtEncseqReader *esr2,
                                        unsigned long pos1,
                                        unsigned long pos2,
                                        unsigned long depth,
                                        unsigned long maxdepth);

/* The following function extracts a twobit encoding at position 
  <pos> with the given <readmode> in the sequence encoded by <encseq>. 
  The <esr> structure refers to a memory area reinitialized in the
  function. The result is stored in <ptbe>. */

void gt_encseq_extract2bitencwithstoppos(GtEndofTwobitencoding *ptbe,
                                         GtEncseqReader *esr,
                                         const GtEncseq *encseq,
                                         GtReadmode readmode,
                                         unsigned long pos);

/* The following function compares the two bit encodings <ptbe1> and <ptbe2>
  and stores the result of the comparison in <commonunits>. The comparison is
  done in forward direction iff <fwd> is true. 
  The comparison is done for the complemented characters
  iff <complement> is true. */
int gt_encseq_compare_pairof_twobitencodings(bool fwd,
                                             bool complement,
                                             GtCommonunits *commonunits,
                                             const GtEndofTwobitencoding *ptbe1,
                                             const GtEndofTwobitencoding 
                                                                      *ptbe2);

/* Return true if and only if the substring of length <len> starting
  at position <startpos> in <encseq> contains a special character.
  <esr> refer to a memory areas for storeing a
  GtEncseqReader. <moveforward> is true if and only if the
  scanning is done in forward direction. */
bool gt_encseq_contains_special(const GtEncseq *encseq,
                                GtReadmode readmode,
                                GtEncseqReader *esr,
                                unsigned long startpos,
                                unsigned long len);

/* Returns the sequence number from the given <position> for an array of of
  SEPARATOR positions <recordseps>.  */
unsigned long gt_encseq_sep2seqnum(const unsigned long *recordseps,
                                   unsigned long numofrecords,
                                   unsigned long totalwidth,
                                   unsigned long position);

/* Returns the number of times that <cc> occurs in the sequences in <encseq>. */
unsigned long gt_encseq_charcount(const GtEncseq *encseq,
                                  GtUchar cc);

/* Prints information about <encseq> via <logger>. If <withfilenames> is set,
  then the filenames of the original sequence files are also printed. */
void gt_encseq_show_features(const GtEncseq *encseq,
                             GtLogger *logger,
                             bool withfilenames);

/* The following generalizes the previous in that the comparison
  of the sequences starts at offset <depth>. */

/* Return the number of positions in <encseq> containing special characters */
unsigned long gt_encseq_specialcharacters(const GtEncseq *encseq);

/* Return the number of ranges of consecutive runs of special characters
  where the length of each range is limited by UCHAR_MAX, USHORT_MAX, and
  UINT32_MAX, depending on whether the GT_ACCESS_TYPE_UCHARTABLES,
  GT_ACCESS_TYPE_USHORTTABLES, GT_ACCESS_TYPE_UINT32TABLES are used */
unsigned long gt_encseq_specialranges(const GtEncseq *encseq);

/* Return the number of ranges of consecutive runs of special characters */
unsigned long gt_encseq_realspecialranges(const GtEncseq *encseq);

/* Return the length of the longest prefix of <encseq> consisting of
  special characters only. */
unsigned long gt_encseq_lengthofspecialprefix(const GtEncseq *encseq);

/* Return the length of the longest suffix of <encseq> consisting of
  special characters only. */
unsigned long gt_encseq_lengthofspecialsuffix(const GtEncseq *encseq);

/* Sets <specialcharinfo> to point to a <GtSpecialcharinfo> for the index
  files specified by <indexname>, even if the encoded sequence is not mapped.
  Returns 0 on success, -1 otherwise. */
int gt_specialcharinfo_read(GtSpecialcharinfo *specialcharinfo,
                            const char *indexname, GtError *err);

/* Returns the encoded representation of the character at position <pos> of
  <encseq> read in the direction as indicated by <readmode>.
  The function only works for sequence representations based on the two bit
  encoding and for the case that encodesequence[pos] does not contain a
  special character. */
GtUchar gt_encseq_extract_encoded_char(const GtEncseq *encseq,
                                       unsigned long pos,
                                       GtReadmode readmode);

/* Sets the sequence input type for <ee> to be pre-encoded. Only for internal
   use. */
void  gt_encseq_encoder_set_input_preencoded(GtEncseqEncoder *ee);

/* The following function shows the encoded sequence at position <startpos>.
   The output goes to the file pointer <fp>. The parameters <fwd> and
   <complement> define whether the sequence is read in forward direction or
   the complement of the sequence is shown. */
void gt_encseq_showatstartpos(FILE *fp,
                              bool fwd,
                              bool complement,
                              const GtEncseq *encseq,
                              unsigned long startpos);

/* Writes the result of gt_encseq_effective_filelength() into the
   location pointed to by <result>. Needed for the Ruby bindings. */
void gt_encseq_effective_filelength_ptr(const GtEncseq *encseq,
                                        uint64_t *result,
                                        unsigned long filenum);

/* The following function shows the encoded sequence at position <startpos> up
   to the first <depth> characters or 30 positions, whichever is the minimum.
   If <depth> is 0, then the entire suffix is shown.
   The output goes to the file pointer <fp>. The parameter readmode and
   define whether the mode in which the sequence is read. */
void gt_encseq_showatstartposwithdepth(FILE *fp,
                                       const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       unsigned long start,
                                       unsigned long depth);

/* Determines the size of the final representation of the encoded sequence
   with the given characteristics, given the access type <sat>. */
uint64_t gt_encseq_determine_size(GtEncseqAccessType sat,
                                  unsigned long totallength,
                                  unsigned long numofdbfiles,
                                  unsigned long lengthofdbfilenames,
                                  unsigned long specialranges,
                                  unsigned int numofchars,
                                  unsigned int bitspersymbol);
/*
  The following function returns the size of the encoded sequence in bytes. */
unsigned long gt_encseq_sizeofrep(const GtEncseq *encseq);

/* The following functions are for testing */

int gt_encseq_builder_unit_test(GtError *err);

/* The following function should only be used for test purposes, because it
  is not efficient. It compares the two suffixes
  at position <start1> and <start2> in <encseq>.  <esr1> and <esr2> refer
  to memory areas for storeing a GtEncseqReader. If <maxdepth>
  is 0, then the entire suffixes are compared. If <maxdepth> is larger than
  0, then only the suffixes up to length <maxdepth> are compared.
  The length of the longest common prefix is stored in <maxlcp>.
  <specialsareequal> specifies if special symbols are considered equal
  during pairwise character comparisons. <specialsareequalatdepth0> specifies
  if special symbols occurring as first symbols of the suffixes
  are considered equal  during pairwise character comparisons.
  The return value is -1, 0 or 1 depending on whether the sequence beginning at
  position <start1> is smaller than, equal to, or larger than the sequence
  beginning at position <start2>. */
int gt_encseq_check_comparetwosuffixes(const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       unsigned long *maxlcp,
                                       bool specialsareequal,
                                       bool specialsareequalatdepth0,
                                       unsigned long maxdepth,
                                       unsigned long start1,
                                       unsigned long start2,
                                       GtEncseqReader *esr1,
                                       GtEncseqReader *esr2);

/* Similar to <gt_error_check()>, this function exits with an error message
  if <encseq> returns inconsistent descriptions when compared to the
  destab. */
void gt_encseq_check_descriptions(const GtEncseq *encseq);

/* Similar to <gt_error_check()>, this function exits with an error message
  if <encseq> returns inconsistent marked positions. */
void gt_encseq_check_markpos(const GtEncseq *encseq);

/* The following function checks the iterators delivering the ranges
  of special characters in an encoded sequence <encseq> in forward and
  in reverse directions. */
int gt_encseq_check_specialranges(const GtEncseq *encseq);

/* The following checks if the encoded sequence <encseq> for consistency.
  It does so by scanning the files given in <filenametab> and comparing
  the extracted symbols to those obtained by directly reading the files.
  Additional, <scantrials> many trials are performed each reading character
  by character starting at some random position and scanning until
  the end of the sequence, while comparing the extraced character
  to the characters extracted by random access. Finally <multicharcmptrials>
  trials are performed each checking the validity of a multicharacter
  extraction.  */
int gt_encseq_check_consistency(const GtEncseq *encseq,
                                const GtStrArray *filenametab,
                                GtReadmode readmode,
                                unsigned long scantrials,
                                unsigned long multicharcmptrials,
                                bool withseqnumcheck,
                                GtError *err);

/* the following function returns a count value */

unsigned long countgt_encseq_compare_viatwobitencoding_get(void);

#endif
