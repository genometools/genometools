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

GtEncodedsequenceOptions* gt_encodedsequence_options_new(void);

void gt_encodedsequence_options_set_progress_timer(GtEncodedsequenceOptions *o,
                                                   GtProgressTimer *pt);
GtProgressTimer* gt_encodedsequence_options_get_progress_timer(
                                                   GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_set_indexname(GtEncodedsequenceOptions *o,
                                              GtStr *indexname);
GtStr* gt_encodedsequence_options_get_indexname(GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_set_symbolmap_file(GtEncodedsequenceOptions *o,
                                                   GtStr *smapfile);
GtStr* gt_encodedsequence_options_get_symbolmap_file(
                                                   GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_set_access_type(GtEncodedsequenceOptions *o,
                                                GtStr *str_sat);
GtStr* gt_encodedsequence_options_get_access_type(GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_set_input_sequences(GtEncodedsequenceOptions *o,
                                                    GtStrArray *filenametab);
GtStrArray* gt_encodedsequence_options_get_input_sequences(
                                                   GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_set_input_dna(GtEncodedsequenceOptions *o);
bool gt_encodedsequence_options_get_input_dna(GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_set_input_protein(GtEncodedsequenceOptions *o);
bool gt_encodedsequence_options_get_input_protein(GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_set_input_plain(GtEncodedsequenceOptions *o);
bool gt_encodedsequence_options_get_input_plain(GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_enable_tis_table_usage(
                                                   GtEncodedsequenceOptions *o);
void gt_encodedsequence_options_disable_tis_table_usage(
                                                   GtEncodedsequenceOptions *o);
bool gt_encodedsequence_options_get_tis_table_usage(
                                                   GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_enable_des_table_usage(
                                                   GtEncodedsequenceOptions *o);
void gt_encodedsequence_options_disable_des_table_usage(
                                                   GtEncodedsequenceOptions *o);
bool gt_encodedsequence_options_get_des_table_usage(
                                                   GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_enable_sds_table_usage(
                                                   GtEncodedsequenceOptions *o);
void gt_encodedsequence_options_disable_sds_table_usage(
                                                   GtEncodedsequenceOptions *o);
bool gt_encodedsequence_options_get_sds_table_usage(
                                                   GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_enable_ssp_table_usage(
                                                   GtEncodedsequenceOptions *o);
void gt_encodedsequence_options_disable_ssp_table_usage(
                                                   GtEncodedsequenceOptions *o);
bool gt_encodedsequence_options_get_ssp_table_usage(
                                                   GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_set_logger(GtEncodedsequenceOptions *o,
                                           GtLogger *logger);
GtLogger* gt_encodedsequence_options_get_logger(GtEncodedsequenceOptions *o);

void gt_encodedsequence_options_delete(GtEncodedsequenceOptions *o);
#endif
