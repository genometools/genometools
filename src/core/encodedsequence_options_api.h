/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#ifndef ENCODEDSEQUENCE_OPTIONS_API_H
#define ENCODEDSEQUENCE_OPTIONS_API_H

#include "core/logger.h"
#include "core/progress_timer.h"
#include "core/str_array.h"

/* The <GtEncodedsequenceOptions> class represents an option set defining
   details for the creation of a <GtEncodedsequence> object. */
typedef struct GtEncodedsequenceOptions GtEncodedsequenceOptions;

/* Creates a new <GtEncodedsequenceOptions> object. */
GtEncodedsequenceOptions* gt_encodedsequence_options_new(void);

/* Sets the progress timer to use during the encoded sequence construction. */
void gt_encodedsequence_options_set_progress_timer(GtEncodedsequenceOptions *o,
                                                   GtProgressTimer *pt);
/* Gets the progress timer registered in <o>. */
GtProgressTimer* gt_encodedsequence_options_get_progress_timer(
                                                   GtEncodedsequenceOptions *o);

/* Sets the index name to use during the encoded sequence construction. */
void gt_encodedsequence_options_set_indexname(GtEncodedsequenceOptions *o,
                                              GtStr *indexname);
/* Gets the index name registered in <o>. */
GtStr* gt_encodedsequence_options_get_indexname(GtEncodedsequenceOptions *o);

/* Sets the symbol map file to use during the encoded sequence construction. */
void gt_encodedsequence_options_set_symbolmap_file(GtEncodedsequenceOptions *o,
                                                   GtStr *smapfile);
/* Gets the symbol map file to registered in <o>. */
GtStr* gt_encodedsequence_options_get_symbolmap_file(
                                                   GtEncodedsequenceOptions *o);

/* Sets the access type  to use during the encoded sequence construction. */
void gt_encodedsequence_options_set_access_type(GtEncodedsequenceOptions *o,
                                                GtStr *str_sat);
/* Gets the access type registered in <o>. */
GtStr* gt_encodedsequence_options_get_access_type(GtEncodedsequenceOptions *o);

/* Sets the string array of input sequence filenames to use during the encoded
   sequence construction. */
void gt_encodedsequence_options_set_input_sequences(GtEncodedsequenceOptions *o,
                                                    GtStrArray *filenametab);
/* Gets the input sequence filenames registered in <o>. */
GtStrArray* gt_encodedsequence_options_get_input_sequences(
                                                   GtEncodedsequenceOptions *o);
/* Declares that the input sequences are DNA sequences. */
void gt_encodedsequence_options_set_input_dna(GtEncodedsequenceOptions *o);
/* Returns true if the input sequences are declared as DNA sequences. */
bool gt_encodedsequence_options_get_input_dna(GtEncodedsequenceOptions *o);

/* Declares that the input sequences are protein sequences. */
void gt_encodedsequence_options_set_input_protein(GtEncodedsequenceOptions *o);
/* Returns true if the input sequences are declared as protein sequences. */
bool gt_encodedsequence_options_get_input_protein(GtEncodedsequenceOptions *o);

/* Declares that the input sequences are pre-encoded sequences. */
void gt_encodedsequence_options_set_input_plain(GtEncodedsequenceOptions *o);
/* Returns true if the input sequences are declared as pre-encoded sequences. */
bool gt_encodedsequence_options_get_input_plain(GtEncodedsequenceOptions *o);

/* Declares that encoded sequences (.esq files) should be created or mapped. */
void gt_encodedsequence_options_enable_tis_table_usage(
                                                   GtEncodedsequenceOptions *o);
/* Declares that encoded sequences (.esq files) should not be created or
   mapped. */
void gt_encodedsequence_options_disable_tis_table_usage(
                                                   GtEncodedsequenceOptions *o);
/* Returns whether encoded sequences (.esq files) should be created or
   mapped. */
bool gt_encodedsequence_options_get_tis_table_usage(
                                                   GtEncodedsequenceOptions *o);

/* Declares that descriptions tables (.des files) should be created or
   mapped. */
void gt_encodedsequence_options_enable_des_table_usage(
                                                   GtEncodedsequenceOptions *o);
/* Declares that descriptions tables (.des files) should not be created or
   mapped. */
void gt_encodedsequence_options_disable_des_table_usage(
                                                   GtEncodedsequenceOptions *o);
/* Returns whether descriptions tables (.des files) should be created or
   mapped. */
bool gt_encodedsequence_options_get_des_table_usage(
                                                   GtEncodedsequenceOptions *o);

/* Declares that description separator positions (.sds files) should be created
   or mapped. */
void gt_encodedsequence_options_enable_sds_table_usage(
                                                   GtEncodedsequenceOptions *o);
/* Declares that description separator positions (.sds files) should not be
   created or mapped. */
void gt_encodedsequence_options_disable_sds_table_usage(
                                                   GtEncodedsequenceOptions *o);
/* Returns whether description separator positions (.sds files) should be
   created or mapped. */
bool gt_encodedsequence_options_get_sds_table_usage(
                                                   GtEncodedsequenceOptions *o);

/* Declares that sequence separator tables (.ssp files) should be created
   or mapped. */
void gt_encodedsequence_options_enable_ssp_table_usage(
                                                   GtEncodedsequenceOptions *o);
/* Declares that sequence separator tables (.ssp files) should not be created
   or mapped. */
void gt_encodedsequence_options_disable_ssp_table_usage(
                                                   GtEncodedsequenceOptions *o);
/* Returns whether sequence separator tables (.ssp files) should be created
   or mapped. */
bool gt_encodedsequence_options_get_ssp_table_usage(
                                                   GtEncodedsequenceOptions *o);

/* Sets the logger object to use during the encoded sequence construction. */
void gt_encodedsequence_options_set_logger(GtEncodedsequenceOptions *o,
                                           GtLogger *logger);
/* Gets the logger object registered in <o>. */
GtLogger* gt_encodedsequence_options_get_logger(GtEncodedsequenceOptions *o);

/* Deletes <o> and frees all associated space. */
void gt_encodedsequence_options_delete(GtEncodedsequenceOptions *o);
#endif
