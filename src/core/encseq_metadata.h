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

#ifndef ENCSEQ_METADATA_H
#define ENCSEQ_METADATA_H

#include "core/alphabet.h"
#include "core/chardef.h"
#include "core/encseq_access_type.h"
#include "core/error_api.h"

typedef struct GtEncseqMetadata GtEncseqMetadata;

GtEncseqMetadata*    gt_encseq_metadata_new(const char *indexname,
                                            GtError *err);
GtAlphabet*          gt_encseq_metadata_alphabet(GtEncseqMetadata *emd);
unsigned long        gt_encseq_metadata_version(GtEncseqMetadata *emd);
bool                 gt_encseq_metadata_is64bit(GtEncseqMetadata *emd);
unsigned long        gt_encseq_metadata_total_length(GtEncseqMetadata *emd);
unsigned long        gt_encseq_metadata_num_of_sequences(GtEncseqMetadata *emd);
unsigned long        gt_encseq_metadata_num_of_files(GtEncseqMetadata *emd);
unsigned long        gt_encseq_metadata_min_seq_length(GtEncseqMetadata *emd);
unsigned long        gt_encseq_metadata_max_seq_length(GtEncseqMetadata *emd);
unsigned long        gt_encseq_metadata_length_of_filenames(
                                                         GtEncseqMetadata *emd);
bool                 gt_encseq_metadata_has_custom_alphabet(
                                                         GtEncseqMetadata *emd);
unsigned long        gt_encseq_metadata_length_of_alphadef(
                                                         GtEncseqMetadata *emd);
GtEncseqAccessType   gt_encseq_metadata_accesstype(GtEncseqMetadata *emd);
GtSpecialcharinfo    gt_encseq_metadata_specialcharinfo(GtEncseqMetadata *emd);
void                 gt_encseq_metadata_delete(GtEncseqMetadata *emd);

#endif
