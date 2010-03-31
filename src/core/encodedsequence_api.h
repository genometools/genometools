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
#include "core/str.h"
#include "core/str_array.h"
#include "core/symboldef.h"

/* Holds information about one sequence in a <GtEncodedsequence>. */
typedef struct
{
  unsigned long seqstartpos,  /* the position of the first character
                                 in the encoded sequence */
                seqlength;    /* the length of the sequence */
} GtSeqinfo;

/* The GtEncodedsequence class represents a collection of sequences from one
   or more input files in a compressed encoding. */
typedef struct GtEncodedsequence GtEncodedsequence;

/* The GtEncodedsequenceScanstate class represents the current state of a
   sequential scan of a <GtEncodedsequence> region. */
typedef struct GtEncodedsequenceScanstate GtEncodedsequenceScanstate;

/* The file suffix used for encoded sequence files. */
#define GT_ENCSEQFILESUFFIX ".esq"
/* The file suffix used for encoded sequence separator position tables. */
#define GT_SSPTABFILESUFFIX ".ssp"
/* The file suffix used for sequence description tables. */
#define GT_DESTABFILESUFFIX ".des"
/* The file suffix used for sequence description separator position tables. */
#define GT_SDSTABFILESUFFIX ".sds"

/* Returns a new <GtEncodedsequence> created from a set of filenames
   (<filenametab>) and creates the on-disk representation on the fly.
   Returns NULL on error.
   The parameter <sfxprogress> specifies an <GtProgressTimer> facilitating log
   state (can be NULL), while <logger> references a <GtLogger> specifying a log
   target.
   <str_indexname> specified a path prefix to use for the newly generated
   table files, creation of which can be toggled via the
   <out{tis,des,sds,ssp}tab> options. The <is{dna,protein,plain}> options
   specify the kind of sequence format given in the source files.
   <str_smap> is the name of an alphabet mapping file declaring an encoding
   scheme for the input sequence(s). <str_sat> is the string representation of
   a position access type and determines how the encoded sequence is stored
   (e.g. "direct", "ushort", "bit", ...). */
GtEncodedsequence* gt_encodedsequence_new_from_files(
                                                  GtProgressTimer *sfxprogress,
                                                  const GtStr *str_indexname,
                                                  const GtStr *str_smap,
                                                  const GtStr *str_sat,
                                                  GtStrArray *filenametab,
                                                  bool isdna,
                                                  bool isprotein,
                                                  bool isplain,
                                                  bool outtistab,
                                                  bool outdestab,
                                                  bool outsdstab,
                                                  bool outssptab,
                                                  GtLogger *logger,
                                                  GtError *err);

/* Returns a new <GtEncodedsequence> created from a set of preprocessed index
   files. Returns NULL on error. The parameter <logger> is used to pass
   a <GtLogger> specifying a log target.
   <indexname> specified a path prefix to the index files to map. Which tables
   are to be mapped can be toggled via the <with{tis,des,sds,ssp}tab> options.
   TODO: document withrange  */
GtEncodedsequence* gt_encodedsequence_new_from_index(bool withrange,
                                                     const GtStr *indexname,
                                                     bool withtistab,
                                                     bool withdestab,
                                                     bool withsdstab,
                                                     bool withssptab,
                                                     GtLogger *logger,
                                                     GtError *err);

/* Returns a new <GtEncodedsequence> created from a pre-encoded sequence in
   memory, given by two sequence pointers <seq1> and <seq2> of lengths <len1>
   and <len2>, respectively. Returns NULL on error.
   <alpha> is the <GtAlphabet> used to encode the sequence. The parameter
   <logger> is used to pass a <GtLogger> specifying a log target.
   TODO: document withrange */
GtEncodedsequence* gt_encodedsequence_new_from_plain(bool withrange,
                                                     const GtUchar *seq1,
                                                     unsigned long len1,
                                                     const GtUchar *seq2,
                                                     unsigned long len2,
                                                     GtAlphabet *alpha,
                                                     GtLogger *logger);

#undef GT_INLINEDENCSEQ
#ifdef GT_INLINEDENCSEQ
  #include "core/encodedsequence_inline.h"
#else

/* Returns the total number of characters in all sequences of <encseq>,
   not including separators. */
unsigned long      gt_encodedsequence_totallength(
                                               const GtEncodedsequence *encseq);

/* Returns the total number of sequences contained in <encseq>. */
unsigned long      gt_encodedsequence_num_of_sequences(
                                               const GtEncodedsequence *encseq);

/* Returns the encoded representation of the character at position <pos> of
   <encseq> read in the direction as indicated by <readmode>. */
GtUchar            gt_encodedsequence_getencodedchar(
                                                const GtEncodedsequence *encseq,
                                                unsigned long pos,
                                                GtReadmode readmode);

/* Returns the encoded representation of the character at position <pos> of
   <encseq> read in the direction as indicated by <readmode>.
   TODO: How is this different from getencodedchar()? */
GtUchar            gt_encodedsequence_extractencodedchar(
                                                const GtEncodedsequence *encseq,
                                                unsigned long pos,
                                                GtReadmode readmode);

/* TODO: please document me */
GtUchar            gt_encodedsequence_getencodedcharnospecial(
                                                const GtEncodedsequence *encseq,
                                                unsigned long pos,
                                                GtReadmode readmode);

/* Returns the encoded representation of the character at position <pos> of
   <encseq> read in the direction as indicated by <readmode>. This function is
   optimized for sequential access to the sequence (e.g. in a for loop). The
   current state of the sequential scan is given by <esr>. */
GtUchar            gt_encodedsequence_sequentialgetencodedchar(
                                                const GtEncodedsequence *encseq,
                                                GtEncodedsequenceScanstate *esr,
                                                unsigned long pos,
                                                GtReadmode readmode);
#endif

/* Returns the encoded representation of the substring from position <frompos>
   to position <topos> of <encseq>. The result is written to the location
   pointed to by <buffer>, which must be large enough to hold the result. */
void               gt_encodedsequence_extract_substring(
                                                const GtEncodedsequence *encseq,
                                                GtUchar *buffer,
                                                unsigned long frompos,
                                                unsigned long topos);

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
GtAlphabet*  gt_encodedsequence_alphabet(const GtEncodedsequence *encseq);

/* Returns a <GtStrArray> of the names of the original sequence files
   contained in <encseq>. */
const GtStrArray*  gt_encodedsequence_filenames(
                                               const GtEncodedsequence *encseq);

/* Deletes <encseq> and frees all associated space. */
void               gt_encodedsequence_delete(GtEncodedsequence *encseq);

/* TODO: please document me */
GtEncodedsequenceScanstate* gt_encodedsequence_scanstate_new(void);

/* TODO: please document me */
void                        gt_encodedsequence_scanstate_init(
                                                GtEncodedsequenceScanstate *esr,
                                                const GtEncodedsequence *encseq,
                                                GtReadmode readmode,
                                                unsigned long startpos);

/* TODO: please document me */
void                        gt_encodedsequence_scanstate_initgeneric(
                                                GtEncodedsequenceScanstate *esr,
                                                const GtEncodedsequence *encseq,
                                                bool moveforward,
                                                unsigned long startpos);

/* TODO: please document me */
void                        gt_encodedsequence_scanstate_delete(
                                               GtEncodedsequenceScanstate *esr);

#endif
