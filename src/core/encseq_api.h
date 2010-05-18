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

#ifndef ENCSEQ_API_H
#define ENCSEQ_API_H

#include "core/alphabet_api.h"
#include "core/logger_api.h"
#include "core/progress_timer_api.h"
#include "core/readmode_api.h"
#include "core/str_api.h"
#include "core/str_array_api.h"

/* The <GtEncseq> class represents a  concatenated collection of sequences from
   one or more input files in a compressed encoding. */
typedef struct GtEncseq GtEncseq;

/* The <GtEncseqEncoder> class creates objects encapsulating a parameter
   set for conversion from sequence files into encoded sequence files on
   secondary storage. */
typedef struct GtEncseqEncoder GtEncseqEncoder;

/* The <GtEncseqLoader> class creates <GtEncseq> objects by mapping index files
   from secondary storage into memory. */
typedef struct GtEncseqLoader GtEncseqLoader;

/* The <GtEncseqBuilder> class creates <GtEncseq> objects by constructing
   uncompressed, encoded string copies in memory. */
typedef struct GtEncseqBuilder GtEncseqBuilder;

/* The <GtEncseqReader> class represents the current state of a
   sequential scan of a <GtEncseq> region as an iterator. */
typedef struct GtEncseqReader GtEncseqReader;

/* The file suffix used for encoded sequence files. */
#define GT_ENCSEQFILESUFFIX ".esq"
/* The file suffix used for encoded sequence separator position tables. */
#define GT_SSPTABFILESUFFIX ".ssp"
/* The file suffix used for sequence description tables. */
#define GT_DESTABFILESUFFIX ".des"
/* The file suffix used for sequence description separator position tables. */
#define GT_SDSTABFILESUFFIX ".sds"

#undef GT_INLINEDENCSEQ
#ifdef GT_INLINEDENCSEQ
  #include "core/encseq_inline.h"
#else
/* Returns the total number of characters in all sequences of <encseq>,
   not including separators. */
unsigned long     gt_encseq_total_length(const GtEncseq *encseq);
/* Returns the total number of sequences contained in <encseq>. */
unsigned long     gt_encseq_num_of_sequences(const GtEncseq *encseq);
/* Returns the encoded representation of the character at position <pos> of
   <encseq> read in the direction as indicated by <readmode>. */
GtUchar           gt_encseq_get_encoded_char(const GtEncseq *encseq,
                                             unsigned long pos,
                                             GtReadmode readmode);
#endif

/* Increases the reference count of <encseq>. */
GtEncseq*         gt_encseq_ref(GtEncseq *encseq);
/* Returns a new <GtEncseqReader> for <encseq>, starting from position
   <startpos>. If <moveforward> is true, the iterator will move forward in the
   sequence, if false, it will move backwards. */
GtEncseqReader*   gt_encseq_create_reader_with_direction(const GtEncseq *encseq,
                                                        bool moveforward,
                                                        unsigned long startpos);
/* Returns a new <GtEncseqReader> for <encseq>, starting from position
   <startpos>. Also supports reading the sequence from the reverse and
   delivering (reverse) complement characters on DNA alphabets using the
   <readmode> option. Please make sure that the <GT_READMODE_COMPL> and
   <GT_READMODE_REVCOMPL> readmodes are only used on DNA alphabets. */
GtEncseqReader*   gt_encseq_create_reader_with_readmode(const GtEncseq *encseq,
                                                        GtReadmode readmode,
                                                        unsigned long startpos);
/* Returns the encoded representation of the substring from position <frompos>
   to position <topos> of <encseq>. The result is written to the location
   pointed to by <buffer>, which must be large enough to hold the result. */
void              gt_encseq_extract_substring(const GtEncseq *encseq,
                                              GtUchar *buffer,
                                              unsigned long frompos,
                                              unsigned long topos);
/* Returns the decoded version of the substring from position <frompos>
   to position <topos> of <encseq>. The result is written to the location
   pointed to by <buffer>, which must be large enough to hold the result. */
void              gt_encseq_extract_decoded(const GtEncseq *encseq,
                                            char *buffer,
                                            unsigned long frompos,
                                            unsigned long topos);
/* Returns the length of the <seqnum>-th sequence in the <encseq>. */
unsigned long     gt_encseq_seqlength(const GtEncseq *encseq,
                                      unsigned long seqnum);
/* Returns the start position of the <seqnum>-th sequence in the <encseq>. */
unsigned long     gt_encseq_seqstartpos(const GtEncseq *encseq,
                                        unsigned long seqnum);
/* Returns a pointer to the description of the <seqnum>-th sequence in the
   <encseq>. The length of the returned string is written to the
   location pointed at by <desclen>. */
const char*       gt_encseq_description(const GtEncseq *encseq,
                                        unsigned long *desclen,
                                        unsigned long seqnum);
/* Returns the <GtAlphabet> associated with <encseq>. */
GtAlphabet*       gt_encseq_alphabet(const GtEncseq *encseq);
/* Returns a <GtStrArray> of the names of the original sequence files
   contained in <encseq>. */
const GtStrArray* gt_encseq_filenames(const GtEncseq *encseq);
/* Deletes <encseq> and frees all associated space. */
void              gt_encseq_delete(GtEncseq *encseq);

/* Reinitializes the given <esr> with the values as described in
   <gt_encseq_create_reader_with_readmode()>. */
void            gt_encseq_reader_reinit_with_readmode(GtEncseqReader *esr,
                                                      const GtEncseq *encseq,
                                                      GtReadmode readmode,
                                                      unsigned long startpos);
/* Reinitializes the given <esr> with the values as described in
   <gt_encseq_create_reader_with_direction()>. */
void            gt_encseq_reader_reinit_with_direction(GtEncseqReader *esr,
                                                       const GtEncseq *encseq,
                                                       bool moveforward,
                                                       unsigned long startpos);
/* Returns the next character from current position of <esr>, advancing the
   iterator by one position. */
GtUchar         gt_encseq_reader_next_encoded_char(GtEncseqReader *esr);
/* Deletes <esr>, freeing all associated space. */
void            gt_encseq_reader_delete(GtEncseqReader *esr);

/* Creates a new <GtEncseqEncoder>. */
GtEncseqEncoder* gt_encseq_encoder_new(void);
/* Sets <pt> to be the progress timer for <ee>. Default is NULL (no progress
   reporting). */
void             gt_encseq_encoder_set_progresstimer(GtEncseqEncoder *ee,
                                                     GtProgressTimer *pt);
/* Sets the representation of <ee> to <sat> which must be one of 'direct',
   'bytecompress', 'bit', 'uchar', 'ushort' or 'uint32'. Returns 0 on success,
   and a negative value on error (<err> is set accordingly). */
int              gt_encseq_encoder_use_representation(GtEncseqEncoder *ee,
                                                      GtStr *sat,
                                                      GtError *err);
/* Sets the symbol map file to use in <ee> to <smap> which must a valid
   alphabet description file. Returns 0 on success, and a negative value on
   error (<err> is set accordingly). Default is NULL (no alphabet
   transformation). */
int              gt_encseq_encoder_use_symbolmap_file(GtEncseqEncoder *ee,
                                                      GtStr *smap,
                                                      GtError *err);
/* Sets the logger to use by <ee> during encoding to <l>. Default is NULL (no
   logging). */
void             gt_encseq_encoder_set_logger(GtEncseqEncoder *ee,
                                              GtLogger *l);
/* Enables support for retrieving descriptions from the encoded sequence
   encoded by <ee>. That is, the .des and .sds tables are created.
   This is a prerequisite for being able to activate description support in
   <gt_encseq_loader_require_description_support()>. Activated by default. */
void             gt_encseq_encoder_enable_description_support(
                                                           GtEncseqEncoder *ee);
/* Disables support for retrieving descriptions from the encoded sequence
   encoded by <ee>. That is, the .des and .sds tables are not created.
   Encoded sequences created without this support will not be able to be
   loaded via a <GtEncseqLoader> with
   <gt_encseq_loader_require_description_support()> enabled. */
void             gt_encseq_encoder_disable_description_support(
                                                           GtEncseqEncoder *ee);
/* Enables support for random access to multiple sequences in the encoded
   sequence encoded by <ee>. That is, the .ssp table is created.
   This is a prerequisite for being able to activate description support in
   <gt_encseq_loader_require_multiseq_support()>. Activated by default. */
void             gt_encseq_encoder_enable_multiseq_support(GtEncseqEncoder *ee);
/* Disables support for random access to multiple sequences in the encoded
   sequence encoded by <ee>. That is, the .ssp table is not created.
   Encoded sequences created without this support will not be able to be
   loaded via a <GtEncseqLoader> with
   <gt_encseq_loader_require_multiseq_support()> enabled. */
void             gt_encseq_encoder_disable_multiseq_support(
                                                           GtEncseqEncoder *ee);
/* Enables creation of the .esq table containing the encoded sequence itself.
   Naturally, enabled by default. */
void             gt_encseq_encoder_create_tis_tab(GtEncseqEncoder *ee);
/* Disables creation of the .esq table. */
void             gt_encseq_encoder_do_not_create_tis_tab(GtEncseqEncoder *ee);
/* Enables creation of the .des table containing sequence descriptions.
   Enabled by default. */
void             gt_encseq_encoder_create_des_tab(GtEncseqEncoder *ee);
/* Disables creation of the .des table. */
void             gt_encseq_encoder_do_not_create_des_tab(GtEncseqEncoder *ee);
/* Enables creation of the .ssp table containing indexes for multiple sequences.
   Enabled by default. */
void             gt_encseq_encoder_create_ssp_tab(GtEncseqEncoder *ee);
/* Disables creation of the .ssp table. */
void             gt_encseq_encoder_do_not_create_ssp_tab(GtEncseqEncoder *ee);
/* Enables creation of the .sds table containing indexes for sequence
   descriptions. Enabled by default. */
void             gt_encseq_encoder_create_sds_tab(GtEncseqEncoder *ee);
/* Disables creation of the .sds table. */
void             gt_encseq_encoder_do_not_create_sds_tab(GtEncseqEncoder *ee);
/* Sets the sequence input type for <ee> to DNA. */
void             gt_encseq_encoder_set_input_dna(GtEncseqEncoder *ee);
/* Sets the sequence input type for <ee> to protein/amino acids. */
void             gt_encseq_encoder_set_input_protein(GtEncseqEncoder *ee);
/* Encodes the sequence files given in <seqfiles> using the settings in <ee>
   and <indexname> as the prefix for the index tables. Returns 0 on success, or
   a negative value on error (<err> is set accordingly). */
int              gt_encseq_encoder_encode(GtEncseqEncoder *ee,
                                          GtStrArray *seqfiles,
                                          const char *indexname,
                                          GtError *err);
/* Deletes <ee>. */
void             gt_encseq_encoder_delete(GtEncseqEncoder *ee);

/* Creates a new <GtEncseqLoader>. */
GtEncseqLoader*  gt_encseq_loader_new(void);
/* Enables support for retrieving descriptions from the encoded sequence
   to be loaded by <el>. That is, the .des and .sds tables must be present.
   For example, these tables are created by having enabled the
   <gt_encseq_encoder_enable_description_support()> option when encoding.
   Activated by default. */
void             gt_encseq_loader_require_description_support(
                                                            GtEncseqLoader *el);
/* Disables support for retrieving descriptions from the encoded sequence
   to be loaded by <el>. That is, the .des and .sds tables need not be present.
   However, disabling this support will result in an error when trying to call
   the method <gt_encseq_description()> on the <GtEncseq>
   object created by <el>. */
void             gt_encseq_loader_drop_description_support(GtEncseqLoader *el);
/* Enables support for random access to multiple sequences in the encoded
   sequence to be loaded by <el>. That is, the .ssp table must be present.
   For example, this table is created by having enabled the
   <gt_encseq_encoder_enable_multiseq_support()> option when encoding.
   Activated by default. */
void             gt_encseq_loader_require_multiseq_support(GtEncseqLoader *el);
/* Disables support for random access to multiple sequences in the encoded
   sequence to be loaded by <el>. That is, the .ssp table needs not be present.
   However, disabling this support will result in an error when trying to call
   the method <gt_encseq_seqinfo()> on the <GtEncseq>
   object created by <el>. */
void             gt_encseq_loader_drop_multiseq_support(GtEncseqLoader *el);
/* Requires presence of the .esq table containing the encoded sequence itself.
   Naturally, enabled by default. */
void             gt_encseq_loader_require_tis_tab(GtEncseqLoader *el);
/* Disables requirement of the .esq table for loading a <GtEncseq>
   using <el>. */
void             gt_encseq_loader_do_not_require_tis_tab(GtEncseqLoader *el);
/* Requires presence of the .des table containing sequence descriptions.
   Enabled by default. */
void             gt_encseq_loader_require_des_tab(GtEncseqLoader *el);
/* Disables requirement of the .des table for loading a <GtEncseq>
   using <el>. */
void             gt_encseq_loader_do_not_require_des_tab(GtEncseqLoader *el);
/* Requires presence of the .ssp table containing indexes for multiple
   sequences. Enabled by default. */
void             gt_encseq_loader_require_ssp_tab(GtEncseqLoader *el);
/* Disables requirement of the .ssp table for loading a <GtEncseq>
   using <el>. */
void             gt_encseq_loader_do_not_require_ssp_tab(GtEncseqLoader *el);
/* Requires presence of the .sds table containing indexes for sequence
   descriptions. Enabled by default. */
void             gt_encseq_loader_require_sds_tab(GtEncseqLoader *el);
/* Disables requirement of the .sds table for loading a <GtEncseq>
   using <el>. */
void             gt_encseq_loader_do_not_require_sds_tab(GtEncseqLoader *el);
/* Enables support for faster access to special ranges in the
   <GtEncseq> loaded by <el>. Enabled by default. */
void             gt_encseq_loader_enable_range_iterator(GtEncseqLoader *el);
/* Disables support for faster access to special ranges in the
   <GtEncseq> loaded by <el>. */
void             gt_encseq_loader_disable_range_iterator(GtEncseqLoader *el);
/* Sets the logger to use by <ee> during encoding to <l>. Default is NULL (no
   logging). */
void             gt_encseq_loader_set_logger(GtEncseqLoader *el, GtLogger *l);
/* Attempts to map the index files as specified by <indexname> using the options
   set in <el> using this interface. Returns a <GtEncseq> instance
   on success, or NULL on error. If an error occurred, <err> is set
   accordingly. */
GtEncseq*        gt_encseq_loader_load(GtEncseqLoader *el,
                                       const char *indexname,
                                       GtError *err);
/* Deletes <el>. */
void             gt_encseq_loader_delete(GtEncseqLoader *el);

/* Creates a new <GtEncseqBuilder> using the alphabet <alpha> as a basis for
   on-the-fly encoding of sequences in memory. */
GtEncseqBuilder* gt_encseq_builder_new(GtAlphabet *alpha);
/* Enables support for faster access to special ranges in the
   <GtEncseq> created by <eb>. Enabled by default. */
void             gt_encseq_builder_enable_range_iterator(GtEncseqBuilder *eb);
/* Disables support for faster access to special ranges in the
   <GtEncseq> created by <eb>. */
void             gt_encseq_builder_disable_range_iterator(GtEncseqBuilder *eb);
/* Enables support for retrieving descriptions from the encoded sequence
   to be built by <eb>. Requires additional memory to hold the descriptions and
   a position index.
   Activated by default. */
void             gt_encseq_builder_enable_description_support(
                                                           GtEncseqBuilder *eb);
/* Disables support for retrieving descriptions from the encoded sequence
   to be built by <eb>. Disabling this support will result in an error when
   trying to call the method <gt_encseq_description()> on the
   <GtEncseq> object created by <eb>. */
void             gt_encseq_builder_disable_description_support(
                                                           GtEncseqBuilder *eb);
/* Enables support for random access to multiple sequences in the encoded
   sequence to be built by <eb>. Requires additional memory for an index of
   starting positions. Activated by default. */
void             gt_encseq_builder_enable_multiseq_support(GtEncseqBuilder *eb);
/* Disables support for random access to multiple sequences in the encoded
   sequence to be built by <eb>. Disabling this support will result in an
   error when trying to call the method <gt_encseq_seqinfo()> on
   the <GtEncseq> object created by <eb> */
void             gt_encseq_builder_disable_multiseq_support(
                                                           GtEncseqBuilder *eb);
/* Enables creation of the .esq table containing the encoded sequence itself.
   Naturally, enabled by default. */
void             gt_encseq_builder_create_tis_tab(GtEncseqBuilder *eb);
/* Disables creation of the .esq table. */
void             gt_encseq_builder_do_not_create_tis_tab(GtEncseqBuilder *eb);
/* Enables creation of the .des table containing sequence descriptions.
   Not enabled by default. */
void             gt_encseq_builder_create_des_tab(GtEncseqBuilder *eb);
/* Disables creation of the .des table. */
void             gt_encseq_builder_do_not_create_des_tab(GtEncseqBuilder *eb);
/* Enables creation of the .ssp table containing indexes for multiple sequences.
   Not enabled by default. */
void             gt_encseq_builder_create_ssp_tab(GtEncseqBuilder *eb);
/* Disables creation of the .ssp table. */
void             gt_encseq_builder_do_not_create_ssp_tab(GtEncseqBuilder *eb);
/* Enables creation of the .sds table containing indexes for sequence
   descriptions. Not enabled by default. */
void             gt_encseq_builder_create_sds_tab(GtEncseqBuilder *eb);
/* Disables creation of the .sds table. */
void             gt_encseq_builder_do_not_create_sds_tab(GtEncseqBuilder *eb);
/* Adds a sequence given as a C string <str> of length <strlen> to the
   encoded sequence to be built by <eb>. Additionally, a description can be
   given (<desc>). If description support is enabled, this must not be NULL.
   A copy will be made during the addition process and the sequence will
   be encoded using the alphabet set at the construction time of <eb>. Thus it
   must only contain symbols compatible with the alphabet. */
void             gt_encseq_builder_add_cstr(GtEncseqBuilder *eb,
                                            const char *str,
                                            unsigned long strlen,
                                            const char *desc);
/* Adds a sequence given as a GtStr <str> to the encoded sequence to be built
   by <eb>. Additionally, a description can be given. If description support
   is enabled, <desc> must not be NULL.
   A copy will be made during the addition process and the sequence will
   be encoded using the alphabet set at the construction time of <eb>. Thus it
   must only contain symbols compatible with the alphabet. */
void             gt_encseq_builder_add_str(GtEncseqBuilder *eb, GtStr *str,
                                           const char *desc);
/* Adds a sequence given as a pre-encoded string  <str> of length <strlen> to
   the encoded sequence to be built by <eb>. <str> must be encoded using the
   alphabet set at the construction time of <eb>.
   Does not take ownership of <str>.
   Additionally, a description <desc> can be given. If description support
   is enabled, this must not be NULL. */
void             gt_encseq_builder_add_encoded(GtEncseqBuilder *eb,
                                               const GtUchar *str,
                                               unsigned long strlen,
                                               const char *desc);
/* Sets the logger to use by <ee> during encoding to <l>. Default is NULL (no
   logging). */
void             gt_encseq_builder_set_logger(GtEncseqBuilder*, GtLogger *l);
/* Creates a new <GtEncseq> from the sequences added to <eb>.
   Returns a <GtEncseq> instance on success, or NULL on error.
   If an error occurred, <err> is set accordingly.
   The state of <eb> is reset to empty after successful creation of a new
   <GtEncseq> (like having called <gt_encseq_builder_reset()>). */
GtEncseq*        gt_encseq_builder_build(GtEncseqBuilder *eb, GtError *err);
/* Clears all added sequences and descriptions, resetting <eb> to a state
   similar to the state immediately after its initial creation.  */
void             gt_encseq_builder_reset(GtEncseqBuilder *eb);
/* Deletes <eb>. */
void             gt_encseq_builder_delete(GtEncseqBuilder *eb);

#endif
