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

#include "core/alphabet_api.h"
#include "core/encodedsequence_options_api.h"
#include "core/logger_api.h"
#include "core/progress_timer_api.h"
#include "core/readmode_api.h"
#include "core/seqinfo.h"
#include "core/str_api.h"
#include "core/str_array_api.h"

/* The <GtEncodedsequence> class represents a collection of sequences from one
   or more input files in a compressed encoding. */
typedef struct GtEncodedsequence GtEncodedsequence;

/* The <GtEncseqEncoder> class creates objects encapsulating a parameter
   set for conversion from sequence files into encoded sequence files. */
typedef struct GtEncseqEncoder GtEncseqEncoder;

/* The <GtEncseqLoader> class creates <GtEncodedsequence> objects
   by mapping index files from disk into memory. */
typedef struct GtEncseqLoader GtEncseqLoader;

/* The <GtEncseqBuilder> class creates <GtEncodedsequence> objects
   by constructing uncompressed, encoded images in memory. */
typedef struct GtEncseqBuilder GtEncseqBuilder;

/* The <GtEncodedsequenceScanstate> class represents the current state of a
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

#undef GT_INLINEDENCSEQ
#ifdef GT_INLINEDENCSEQ
  #include "core/encodedsequence_inline.h"
#else

/* Returns the total number of characters in all sequences of <encseq>,
   not including separators. */
unsigned long      gt_encodedsequence_total_length(
                                               const GtEncodedsequence *encseq);

/* Returns the total number of sequences contained in <encseq>. */
unsigned long      gt_encodedsequence_num_of_sequences(
                                               const GtEncodedsequence *encseq);

/* Returns the encoded representation of the character at position <pos> of
   <encseq> read in the direction as indicated by <readmode>. */
GtUchar            gt_encodedsequence_get_encoded_char(
                                                const GtEncodedsequence *encseq,
                                                unsigned long pos,
                                                GtReadmode readmode);

/* Returns the encoded representation of the character at position <pos> of
   <encseq> read in the direction as indicated by <readmode>. This function is
   optimized for sequential access to the sequence (e.g. in a for loop). The
   current state of the sequential scan is maintained in <esr>. */
GtUchar            gt_encodedsequence_get_encoded_char_sequential(
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
GtAlphabet*        gt_encodedsequence_alphabet(const GtEncodedsequence *encseq);

/* Returns a <GtStrArray> of the names of the original sequence files
   contained in <encseq>. */
const GtStrArray*  gt_encodedsequence_filenames(
                                               const GtEncodedsequence *encseq);

/* Deletes <encseq> and frees all associated space. */
void               gt_encodedsequence_delete(GtEncodedsequence *encseq);

/* Returns a new object encapsulating the current state of a sequential
   <GtEncodedsequence> scan. */
GtEncodedsequenceScanstate* gt_encodedsequence_scanstate_new_empty(void);

/* Returns a new object encapsulating the current state of a sequential
   <GtEncodedsequence> scan, optimising access to <encseq> starting at position
   <startpos> in the direction specified in <readmode>. */
GtEncodedsequenceScanstate* gt_encodedsequence_scanstate_new(
                                                const GtEncodedsequence *encseq,
                                                GtReadmode readmode,
                                                unsigned long startpos);

/* Reinitializes the given <esr> with the values as described in
   <gt_encodedsequence_scanstate_new()>. */
void                        gt_encodedsequence_scanstate_init(
                                                GtEncodedsequenceScanstate *esr,
                                                const GtEncodedsequence *encseq,
                                                GtReadmode readmode,
                                                unsigned long startpos);

/* Deletes <esr>, freeing all associated space. */
void                        gt_encodedsequence_scanstate_delete(
                                               GtEncodedsequenceScanstate *esr);

GtEncseqEncoder* gt_encseq_encoder_new(void);

void             gt_encseq_encoder_set_progresstimer(GtEncseqEncoder *ee,
                                                     GtProgressTimer *pt);
int              gt_encseq_encoder_use_representation(GtEncseqEncoder *ee,
                                                      GtStr *sat,
                                                      GtError *err);
int              gt_encseq_encoder_use_symbolmap_file(GtEncseqEncoder *ee,
                                                      GtStr *smap,
                                                      GtError *err);
void             gt_encseq_encoder_set_logger(GtEncseqEncoder *ee,
                                              GtLogger *l);

void             gt_encseq_encoder_enable_description_support(
                                                           GtEncseqEncoder *ee);
void             gt_encseq_encoder_disable_description_support(
                                                           GtEncseqEncoder *ee);

void             gt_encseq_encoder_enable_multiseq_support(GtEncseqEncoder *ee);
void             gt_encseq_encoder_disable_multiseq_support(
                                                           GtEncseqEncoder *ee);

void             gt_encseq_encoder_create_tis_tab(GtEncseqEncoder *ee);
void             gt_encseq_encoder_do_not_create_tis_tab(GtEncseqEncoder *ee);

void             gt_encseq_encoder_create_des_tab(GtEncseqEncoder *ee);
void             gt_encseq_encoder_do_not_create_des_tab(GtEncseqEncoder *ee);

void             gt_encseq_encoder_create_ssp_tab(GtEncseqEncoder *ee);
void             gt_encseq_encoder_do_not_create_ssp_tab(GtEncseqEncoder *ee);

void             gt_encseq_encoder_create_sds_tab(GtEncseqEncoder *ee);
void             gt_encseq_encoder_do_not_create_sds_tab(GtEncseqEncoder *ee);

void             gt_encseq_encoder_set_input_dna(GtEncseqEncoder *ee);
void             gt_encseq_encoder_set_input_protein(GtEncseqEncoder *ee);
void             gt_encseq_encoder_set_input_preencoded(GtEncseqEncoder *ee);

int              gt_encseq_encoder_encode(GtEncseqEncoder *ee,
                                          GtStrArray *seqfiles,
                                          GtStr *indexname,
                                          GtError *err);

void             gt_encseq_encoder_delete(GtEncseqEncoder *ee);

GtEncseqLoader*       gt_encseq_loader_new(void);

void                  gt_encseq_loader_require_description_support(
                                                            GtEncseqLoader *el);
void                  gt_encseq_loader_drop_description_support(
                                                            GtEncseqLoader *el);

void                  gt_encseq_loader_require_multiseq_support(
                                                            GtEncseqLoader *el);
void                  gt_encseq_loader_drop_multiseq_support(
                                                            GtEncseqLoader *el);

void                  gt_encseq_loader_require_tis_tab(GtEncseqLoader *el);
void                  gt_encseq_loader_do_not_require_tis_tab(
                                                            GtEncseqLoader *el);

void                  gt_encseq_loader_require_des_tab(GtEncseqLoader *el);
void                  gt_encseq_loader_do_not_require_des_tab(
                                                            GtEncseqLoader *el);

void                  gt_encseq_loader_require_ssp_tab(GtEncseqLoader *el);
void                  gt_encseq_loader_do_not_require_ssp_tab(
                                                            GtEncseqLoader *el);

void                  gt_encseq_loader_require_sds_tab(GtEncseqLoader *el);
void                  gt_encseq_loader_do_not_require_sds_tab(
                                                            GtEncseqLoader *el);

void                  gt_encseq_loader_enable_range_iterator(
                                                            GtEncseqLoader *el);
void                  gt_encseq_loader_disable_range_iterator(
                                                            GtEncseqLoader *el);

void                  gt_encseq_loader_set_logger(GtEncseqLoader *el,
                                                  GtLogger *l);

GtEncodedsequence*    gt_encseq_loader_load(GtEncseqLoader *el,
                                            GtStr *indexname,
                                            GtError *err);

void                  gt_encseq_loader_delete(GtEncseqLoader *el);

GtEncseqBuilder*      gt_encseq_builder_new(GtAlphabet *alpha);

void                  gt_encseq_builder_enable_range_iterator(
                                                           GtEncseqBuilder *eb);
void                  gt_encseq_builder_disable_range_iterator(
                                                           GtEncseqBuilder *eb);

void                  gt_encseq_builder_enable_description_support(
                                                           GtEncseqBuilder *eb);
void                  gt_encseq_builder_disable_description_support(
                                                           GtEncseqBuilder *eb);

void                  gt_encseq_builder_enable_multiseq_support(
                                                           GtEncseqBuilder *eb);
void                  gt_encseq_builder_disable_multiseq_support(
                                                           GtEncseqBuilder *eb);

void                  gt_encseq_builder_create_tis_tab(GtEncseqBuilder *eb);
void                  gt_encseq_builder_do_not_create_tis_tab(
                                                           GtEncseqBuilder *eb);

void                  gt_encseq_builder_create_des_tab(GtEncseqBuilder *eb);
void                  gt_encseq_builder_do_not_create_des_tab(
                                                           GtEncseqBuilder *eb);

void                  gt_encseq_builder_create_ssp_tab(GtEncseqBuilder *eb);
void                  gt_encseq_builder_do_not_create_ssp_tab(
                                                           GtEncseqBuilder *eb);

void                  gt_encseq_builder_create_sds_tab(GtEncseqBuilder *eb);
void                  gt_encseq_builder_do_not_create_sds_tab(
                                                           GtEncseqBuilder *eb);

void                  gt_encseq_builder_add_cstr(GtEncseqBuilder *eb,
                                                 const char *str,
                                                 unsigned long strlen,
                                                 const char *desc);

void                  gt_encseq_builder_add_str(GtEncseqBuilder *eb,
                                                GtStr *str,
                                                const char *desc);

void                  gt_encseq_builder_add_encoded(GtEncseqBuilder *eb,
                                                    const GtUchar *str,
                                                    unsigned long strlen,
                                                    const char *desc);
/* default: NULL */
void                  gt_encseq_builder_set_logger(GtEncseqBuilder*,
                                                   GtLogger *l);

GtEncodedsequence*    gt_encseq_builder_build(GtEncseqBuilder *eb,
                                              GtError *err);

void                  gt_encseq_builder_reset(GtEncseqBuilder *eb);

void                  gt_encseq_builder_delete(GtEncseqBuilder *eb);
#endif
