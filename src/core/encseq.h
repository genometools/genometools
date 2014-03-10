/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2007-2011 Center for Bioinformatics, University of Hamburg

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
#include "core/disc_distri_api.h"
#include "core/encseq_api.h"
#include "core/encseq_access_type.h"
#include "core/encseq_options.h"
#include "core/filelengthvalues.h"
#include "core/intbits.h"
#include "core/md5_tab.h"
#include "core/range.h"
#include "core/readmode.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/types_api.h"
#include "core/arraydef.h"

#define GT_ENCSEQ_VERSION  3

#define GT_REVERSEPOS(TOTALLENGTH,POS) \
          ((TOTALLENGTH) - 1 - (POS))

GT_DECLAREARRAYSTRUCT(GtTwobitencoding);

/* Creates a new <GtEncseqEncoder> using the options given in <opts>.
   If no encoder could be created using the given options, <NULL> is
   returned and <err> is set accordingly. */
GtEncseqEncoder* gt_encseq_encoder_new_from_options(GtEncseqOptions *opts,
                                                    GtError *err);

/* Creates a new <GtEncseqLoader> using the options given in <opts>.
   If no loader could be created using the given options, <NULL> is
   returned and <err> is set accordingly. */
GtEncseqLoader* gt_encseq_loader_new_from_options(GtEncseqOptions *opts,
                                                  GtError *err);

/* Do not output the header of the esq-file. This is needed
   for generating testcases for Timo Beller. */

void gt_encseq_encoder_disable_esq_header(GtEncseqEncoder *ee);

/* The following type stores a two bit encoding in <tbe> with information
  about the number of two bit units which do not store a special
  character in <unitsnotspecial>. To allow the comparison of these
  structures, the <position> in the sequence at which the encoding was
  extracted is also stored */
typedef struct
{
  GtTwobitencoding tbe;           /* two bit encoding */
  unsigned int unitsnotspecial;   /* units which are not special */
  GtUword referstartpos;    /* position of suffix to be compared */
} GtEndofTwobitencoding;

/* The following type stores the result of comparing a pair of twobit
  encodings. <common> stores the number of units which are common
  (either from the beginning or from the end. common is in the range 0 to
  GT_UNITSIN2BITENC. <leftspecial> is true
  iff the first twobitencoding contains a special character in the units
  it covers. <rightspecial> is true iff the second twobit encoding contains
  a special character in the units it covers. <finaldepth> stores
  result of the comparison, i.e. the longest common prefix of the
  comparison. */
typedef struct
{
  unsigned int common;
  bool leftspecial,
       rightspecial;
  GtUword finaldepth;
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
                                             GtUword pos,
                                             GtReadmode readmode);

/* The following function extracts from the byte sequence of length <len>
  pointed to by <seq> the sequence of len/4 bytecodes each
  encoding four consecutive characters in one byte. */
void gt_encseq_plainseq2bytecode(GtUchar *bytecode,
                                 const GtUchar *seq,
                                 GtUword len);

/* The following function extracts from an encoded sequence a substring of
  length <len> beginning at position <startindex> and stores the result
  in the byte sequence encoding four consecutive characters in one byte. */
void gt_encseq_sequence2bytecode(GtUchar *dest,
                                 const GtEncseq *encseq,
                                 GtUword startindex,
                                 GtUword len);

/* Returns true is <encseq> has special ranges, false otherwise. */
bool gt_encseq_has_specialranges(const GtEncseq *encseq);

/* Return true if the representation of the <encseq> is based on two
  bit encoding */
bool gt_encseq_bitwise_cmp_ok(const GtEncseq *encseq);

/* Return the integer code of the sequence of length <prefixlength>
  beginning at position <frompos> represented by <encseq>. <esr> is
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
                                       GtUword frompos,
                                       unsigned int prefixlength);

/* The following function compares two substrings in <encseq1> and
  <encseq2> beginning at position <pos1>+<depth> and <pos2>+<depth>, resp.
  <esr1> and <esr2> refer to memory areas for storeing a GtEncseqReader.
  The information about the length of the longest common prefix is stored
  in <commonunits>. <fwd> and <complement> specify if the sequence is
  scanned in forward direction and if the complement of the sequence is to
  be considered. Thereturn value is -1, 0 or 1 depending on whether the
  sequence beginning at position <pos1>+<depth> is smaller than, equal to,
  or larger than the sequence beginning at position <pos2>+<depth>.
  If <madepth> is 0, then the entire suffixes are compared. Otherwise, the
  comparison is restricted to the prefixes of length <maxdepth>. */
int gt_encseq_compare_viatwobitencoding(GtCommonunits *commonunits,
                                        const GtEncseq *encseq1,
                                        const GtEncseq *encseq2,
                                        GtReadmode readmode,
                                        GtEncseqReader *esr1,
                                        GtEncseqReader *esr2,
                                        GtUword pos1,
                                        GtUword pos2,
                                        GtUword depth,
                                        GtUword maxdepth);

const GtTwobitencoding *gt_encseq_twobitencoding_export(const GtEncseq *encseq);

/* The following two _mapoffset functions are for internal use only.
 * They return an offset value by which the sequence (i.e. the twobitencoding,
 * the plainseq or the bitpackarray) and the
 * chardistri are found in the encseq representation, thus allowing to
 * modify in place the representation by mapping the esq file. */
size_t gt_encseq_sequence_mapoffset(const GtEncseq *encseq);
size_t gt_encseq_chardistri_mapoffset(const GtEncseq *encseq);

/* Saves an encoded sequence characterized by the given parameters into the
   index of the given name. Only sequence collections of type eqlength are
   supported. The other parameters must be consistent with this. */
int gt_encseq_equallength_write_twobitencoding_to_file(const char *indexname,
                                     GtUword totallength,
                                     GtUword lengthofsinglesequence,
                                     GtTwobitencoding *twobitencoding,
                                     GtUword numofsequences,
                                     GtUword numoffiles,
                                     const GtFilelengthvalues *filelengthtab,
                                     const GtStrArray *filenametab,
                                     const GtUword *characterdistribution,
                                     GtError *err);

/* Saves an encoded sequence characterized by the given parameters into the
   index of the given name. Only sequence collections of type
   are UCHAR, USHORT, UINT32 are supported. The other parameters must be
   consistent with this. */
int gt_encseq_generic_write_twobitencoding_to_file(const char *indexname,
                                     GtUword totallength,
                                     GtEncseqAccessType sat,
                                     GtUword lengthofsinglesequence,
                                     GtUword minseqlen,
                                     GtUword maxseqlen,
                                     GtUword lengthofspecialprefix,
                                     GtUword lengthofspecialsuffix,
                                     GtUword lengthoflongestnonspecial,
                                     GtTwobitencoding *twobitencoding,
                                     GtUword numofsequences,
                                     GtUword numoffiles,
                                     const GtFilelengthvalues *filelengthtab,
                                     const GtStrArray *filenametab,
                                     const GtUword *characterdistribution,
                                     GtError *err);

/* The following type is used for computing stoppositions when
   the twobitencoding is used */
typedef struct GtViatwobitkeyvalues GtViatwobitkeyvalues;

/* The following is the constructor for the latter type */
GtViatwobitkeyvalues *gt_Viatwobitkeyvalues_new(void);

/* The following reinitializes the latter type */
void gt_Viatwobitkeyvalues_reinit(GtViatwobitkeyvalues *vtk,
                                  const GtEncseq *encseq,
                                  GtReadmode readmode,
                                  GtEncseqReader *esr,
                                  GtUword pos,
                                  GtUword depth,
                                  GtUword maxdepth,
                                  GtUword stoppos);

/* The following is the destructor for the latter type */
void gt_Viatwobitkeyvalues_delete(GtViatwobitkeyvalues *vtk);

/* The following function compare two substrings of the given encseqs
   <encseq1> and <encseq2>
   refered to via the <GtViatwobitkeyvalues>-parameters. The comparison
   only works for sequences represented as a two bit encodings. It
   starts with offset <depth> and stops at offset <maxdepth>, whenever
   <maxdepth> is larger than 0. If <maxdepth> is 0, then the comparison
   stops at the end of the sequence or at the first special character.
   The encoded sequence is given via the parameter <encseq> and it is
   accessed via the given <readmode> which also determines if the longest
   common prefix is computed (whenever <GT_ISDIRREVERSE(readmode)> is false).
*/
int gt_encseq_twobitencoding_strcmp(GtCommonunits *commonunits,
                                    const GtEncseq *encseq1,
                                    const GtEncseq *encseq2,
                                    GtReadmode readmode,
                                    GtUword depth,
                                    GtUword maxdepth,
                                    GtViatwobitkeyvalues *vtk1,
                                    GtViatwobitkeyvalues *vtk2);

bool gt_encseq_has_twobitencoding(const GtEncseq *encseq);

bool gt_encseq_has_twobitencoding_stoppos_support(const GtEncseq *encseq);

GtUword gt_getnexttwobitencodingstoppos(bool fwd, GtEncseqReader *esr);

/* The following function extracts a twobit encoding at position
  <pos> with the given <readmode> in the sequence encoded by <encseq>.
  The <esr> structure refers to a memory area reinitialized in the
  function. The result is stored in <ptbe>. */

GtUword gt_encseq_extract2bitencwithtwobitencodingstoppos(
                                         GtEndofTwobitencoding *ptbe,
                                         GtEncseqReader *esr,
                                         const GtEncseq *encseq,
                                         GtReadmode readmode,
                                         GtUword pos);

/* The following function extracts the twobitencoding beginning at position
   <pos> in the given encoded sequence <encseq> and ending at
   the next stop position wrt. to the readmode. The extraction
   stops after at most <maxdepth> character have been extracted.
   The result is stored in <tbereservoir> which points to an array of
   size <sizeofvector>. */

unsigned int gt_encseq_extract2bitencvector(
                                         GtArrayGtTwobitencoding *tbereservoir,
                                         const GtEncseq *encseq,
                                         GtEncseqReader *esr,
                                         GtReadmode readmode,
                                         GtUword pos,
                                         bool withstoppos,
                                         GtUword stoppos);

/* The following function extracts the twobitencoding for the
   sequence with sequence number <seqnum> beginning at the relative
   position <relpos>. If <maxnofelem> > 0, then at most <maxnofelem> elements
   are extracted. The result is stored in <tbereservoir> which
   points to an array of appropriate size. The number of elements
   is stored at the address <storedvalues> points to.  */

unsigned int gt_encseq_relpos_extract2bitencvector(
                                          GtArrayGtTwobitencoding *tbereservoir,
                                          const GtEncseq *encseq,
                                          GtUword seqnum,
                                          GtUword relpos,
                                          GtUword maxnofelem);

/* The following function compares the two bit encodings <ptbe1> and <ptbe2>
  and stores the result of the comparison in <commonunits>. The comparison is
  done in forward direction iff <fwd> is true.
  The comparison is done for the complemented characters iff <complement>
  is true. */
int gt_encseq_compare_pairof_twobitencodings(bool fwd,
                                           bool complement,
                                           GtCommonunits *commonunits,
                                           const GtEndofTwobitencoding *ptbe1,
                                           const GtEndofTwobitencoding *ptbe2);

/* The following function compares the two bit encodings and returns the
   length of the longest common prefix of the sequence they represent */

unsigned int gt_encseq_lcpofdifferenttwobitencodings(GtTwobitencoding tbe1,
                                                      GtTwobitencoding tbe2);

/* The following function compares the two bit encodings <tbe1> and <tbe2>
  and stores the result of the comparison in <commonunits>. The comparison is
  done in forward direction iff <fwd> is true.
  The comparison is done for the complemented characters iff <complement>
  is true. */

int gt_encseq_compare_pairof_different_twobitencodings(
                                            bool fwd,
                                            bool complement,
                                            GtCommonunits *commonunits,
                                            GtTwobitencoding tbe1,
                                            GtTwobitencoding tbe2);

/* Return true if and only if the substring of length <len> starting
  at position <startpos> in <encseq> contains a special character.
  <esr> refers to a memory area for storing a GtEncseqReader. */
bool gt_encseq_contains_special(const GtEncseq *encseq,
                                GtReadmode readmode,
                                GtEncseqReader *esr,
                                GtUword startpos,
                                GtUword len);

/* Returns the sequence number from the given <position> for an array of
  SEPARATOR positions <recordseps>.  */
GtUword gt_encseq_sep2seqnum(const GtUword *recordseps,
                                   GtUword numofrecords,
                                   GtUword totalwidth,
                                   GtUword position);

/* Returns the number of times that <cc> occurs in the sequences in <encseq>. */
GtUword gt_encseq_charcount(const GtEncseq *encseq,
                                  GtUchar cc);

/* Prints information about <encseq> via <logger>. If <withfilenames> is set,
  then the filenames of the original sequence files are also printed. */
void gt_encseq_show_features(const GtEncseq *encseq,
                             GtLogger *logger,
                             bool withfilenames);

/* The following generalizes the previous in that the comparison
  of the sequences starts at offset <depth>. */

/* Return the number of positions in <encseq> containing special characters */
GtUword gt_encseq_specialcharacters(const GtEncseq *encseq);

/* Return the number of positions in <encseq> containing wildcard */
GtUword gt_encseq_wildcards(const GtEncseq *encseq);

/* Return the number of ranges of consecutive runs of special characters
  where the length of each range is limited by UCHAR_MAX, USHORT_MAX, and
  UINT32_MAX, depending on whether the GT_ACCESS_TYPE_UCHARTABLES,
  GT_ACCESS_TYPE_USHORTTABLES, GT_ACCESS_TYPE_UINT32TABLES are used */
GtUword gt_encseq_specialranges(const GtEncseq *encseq);

GtUword gt_encseq_exceptioncharacters(const GtEncseq *encseq);

GtUword gt_encseq_exceptionranges(const GtEncseq *encseq);

GtUword gt_encseq_max_subalpha_size(const GtEncseq *encseq);

/* Return the number of ranges of consecutive runs of wildcards
  where the length of each range is limited by UCHAR_MAX, USHORT_MAX, and
  UINT32_MAX, depending on whether the GT_ACCESS_TYPE_UCHARTABLES,
  GT_ACCESS_TYPE_USHORTTABLES, GT_ACCESS_TYPE_UINT32TABLES are used */
GtUword gt_encseq_wildcardranges(const GtEncseq *encseq);

/* Return the number of ranges of consecutive runs of special characters */
GtUword gt_encseq_realspecialranges(const GtEncseq *encseq);

/* Return the number of ranges of consecutive runs of wildcards */
GtUword gt_encseq_realwildcardranges(const GtEncseq *encseq);

/* Return the length of the longest consecutive run of non-special characters */
GtUword gt_encseq_lengthoflongestnonspecial(const GtEncseq *encseq);

/* Return the length of the longest prefix of <encseq> consisting of
  special characters only. */
GtUword gt_encseq_lengthofspecialprefix(const GtEncseq *encseq);

/* Return the length of the longest prefix of <encseq> consisting of
  wildcards only. */
GtUword gt_encseq_lengthofwildcardprefix(const GtEncseq *encseq);

/* Return the length of the longest suffix of <encseq> consisting of
  special characters only. */
GtUword gt_encseq_lengthofspecialsuffix(const GtEncseq *encseq);

/* Return the length of the longest suffix of <encseq> consisting of
  wildcards only. */
GtUword gt_encseq_lengthofwildcardsuffix(const GtEncseq *encseq);

/* Returns number of characters in the alphabet which is part of the
   <encseq>. The number does not include the wildcards. */
unsigned int  gt_encseq_alphabetnumofchars(const GtEncseq *encseq);

/* Returns an array of the characters in the alphabet for the encoded
   sequence <encseq>. */
const GtUchar *gt_encseq_alphabetcharacters(const GtEncseq *encseq);

/* Sets <specialcharinfo> to point to a <GtSpecialcharinfo> for the index
  files specified by <indexname>, even if the encoded sequence is not mapped.
  Returns 0 on success, -1 otherwise. */
int gt_specialcharinfo_read(GtSpecialcharinfo *specialcharinfo,
                            const char *indexname, GtError *err);

/* Sets the sequence input type for <ee> to be pre-encoded. Only for internal
   use. */
void  gt_encseq_encoder_set_input_preencoded(GtEncseqEncoder *ee);

/* Returns <true> if the input sequence has been defined as being pre-encoded.
 */
bool gt_encseq_encoder_is_input_preencoded(GtEncseqEncoder *ee);

/* The following function shows the encoded sequence at position <startpos>.
   The output goes to the file pointer <fp>. The parameters <fwd> and
   <complement> define whether the sequence is read in forward direction or
   the complement of the sequence is shown. */
void gt_encseq_showatstartpos(FILE *fp,
                              bool fwd,
                              bool complement,
                              const GtEncseq *encseq,
                              GtUword startpos);

/* The following function shows the encoded sequence at position <startpos> up
   to the first <depth> characters or 30 positions, whichever is the minimum.
   If <depth> is 0, then the entire suffix is shown.
   The output goes to the file pointer <fp>. The parameter readmode and
   define whether the mode in which the sequence is read. */
void gt_encseq_showatstartposwithdepth(FILE *fp,
                                       const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       GtUword start,
                                       GtUword depth);

/* Determines the size of the final representation of the encoded sequence
   with the given characteristics, given the access type <sat>. */
uint64_t gt_encseq_determine_size(GtEncseqAccessType sat,
                                  GtUword totallength,
                                  GtUword numofsequences,
                                  GtUword numofdbfiles,
                                  GtUword lengthofdbfilenames,
                                  GtUword wildcardranges,
                                  unsigned int numofchars,
                                  unsigned int bitspersymbol,
                                  GtUword lengthofalphadef);

/* The following function returns the size of the encoded sequence in bytes. */
GtUword gt_encseq_sizeofrep(const GtEncseq *encseq);

/* The following function returns the size of the GtEncseq-structure
   in bytes. */
GtUword gt_encseq_sizeofstructure(void);

/* The following function delivers the accesstype of a given encoded
   sequence. */
GtEncseqAccessType gt_encseq_accesstype_get(const GtEncseq *encseq);

/* The following function delivers the encseq->equallength.valueunsignedlong
   if encseq->equallength.defined is true */
GtUword gt_encseq_equallength(const GtEncseq *encseq);

/* Returns a new <GtMD5Tab> allowing to match MD5 sums to sequence numbers for a
   given <encseq>. */
GtMD5Tab*  gt_encseq_get_md5_tab(const GtEncseq *encseq, GtError *err);

/* for a given array of at least one separator positions */
int gt_encseq_seppos2ssptab(const char *indexname,
                            GtUword totallength,
                            GtUword numofdbsequences,
                            const GtUword *seppostab,
                            GtError *err);

/* The following functions are for testing */

#ifndef NDEBUG
void gt_GtSpecialcharinfo_check(const GtSpecialcharinfo *specialcharinfo,
                                GtUword numofseparatorpositions);
#endif

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
  if special symbols occurring as first symbols of the suffixes are considered
  equal during pairwise character comparisons.
  The return value is -1, 0 or 1 depending on whether the sequence beginning at
  position <start1> is smaller than, equal to, or larger than the sequence
  beginning at position <start2>. */
int gt_encseq_check_comparetwosuffixes(const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       GtUword *maxlcp,
                                       bool specialsareequal,
                                       bool specialsareequalatdepth0,
                                       GtUword maxdepth,
                                       GtUword start1,
                                       GtUword start2,
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
void gt_encseq_check_specialranges(const GtEncseq *encseq);

/* Checks whether the information given by gt_encseq_seqstartpos() agrees
   with the actual positions of separators in the encoded sequence.
   A logger allows to report the state of the computation. */
void gt_encseq_check_startpositions(const GtEncseq *encseq,GtLogger *logger);

/* Checks whether the minima/maxima given by gt_encseq_(min,max)_seq_length()
   agree with the sequence lengths in the encoded sequence. */
int gt_encseq_check_minmax(const GtEncseq *encseq, GtError *err);

/* The following checks if the encoded sequence <encseq> for consistency.
  It does so by scanning the files given in <filenametab> and comparing
  the extracted symbols to those obtained by directly reading the files.
  Additional, <scantrials> many trials are performed each reading character
  by character starting at some random position and scanning until
  the end of the sequence, while comparing the extracted character
  to the characters extracted by random access. Finally <multicharcmptrials>
  trials are performed each checking the validity of a multicharacter
  extraction.  */
int gt_encseq_check_consistency(const GtEncseq *encseq,
                                const GtStrArray *filenametab,
                                GtReadmode readmode,
                                GtUword scantrials,
                                GtUword multicharcmptrials,
                                bool withseqnumcheck,
                                bool withcheckunit,
                                GtLogger *logger,
                                GtError *err);

/* The following checks if the function to write an own encoded
   sequence to file works correctly. Requires that the index given by
   <indexname> has an access type of GT_ACCESS_TYPE_EQUALLENGTH. */
int gt_encseq_check_external_twobitencoding_to_file(const char *indexname,
                                                    GtError *err);

uint64_t gt_encseq_pairbitsum(const GtEncseq *encseq);

/* the following function checks if a sequence number and relative
   position is consistent with the given position. It only works for
   sequences of equal length and should be extended to sequences
   of variable length. */

void gt_encseq_relpos_seqnum_check(const char *filename,int line,
                                   const GtEncseq *encseq,GtUword relpos,
                                   GtUword seqnum,GtUword position);

/* for a position outside the range from 0 to totallength -1 deliver a
   unique integer */
#define GT_UNIQUEINT(POS)        ((GtUword) ((POS) + GT_COMPAREOFFSET))
#define GT_ISUNIQUEINT(POS)      ((POS) >= GT_COMPAREOFFSET)

/* Reverse the range with respect to the given total length */
void gt_range_reverse(GtUword totallength,GtRange *range);

/* Return the length of the longest description in <encseq>.
   Requires that the description support is enabled in <encseq>.  */
GtUword gt_encseq_max_desc_length(const GtEncseq *encseq);

#endif
