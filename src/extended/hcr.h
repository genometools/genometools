/*
  Copyright (c) 2011 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef HCR_H
#define HCR_H

#include "core/alphabet_api.h"
#include "core/error.h"
#include "core/str_api.h"
#include "core/timer_api.h"
#include "extended/huffcode.h"
#include "extended/sampling.h"

#define HCRFILESUFFIX ".hcr"
#define HCRFILEDECODEDSUFFIX ".fastq"

/* The classes <GtHcrEncoder> and <GtHcrDecoder> can be used to compress or
   decompress fastq-files containing short reads */
typedef struct GtHcrEncoder GtHcrEncoder;
typedef struct GtHcrDecoder GtHcrDecoder;

/* <GtQualRange> is a struct to define a range of fastq qualities that will be
   kept while compressing fastq short reads */
typedef struct GtQualRange GtQualRange;

struct GtQualRange {
  unsigned start,
           end;
};

/* Returns a new GtHcrEncoder object. If <descs> is true, description lines
   will be encoded. <qrange> denotes the range of quality values. All quality
   values smaller or equal to the lower bound will be converted to the lower
   bound. All quality values equal or larger than the upper bound will be
   converted to the upper bound.
   Reads have to be of constant length, this might be changed in future updates
*/
GtHcrEncoder* gt_hcr_encoder_new(GtStrArray *files, GtAlphabet *alpha,
                                 bool descs, GtQualRange qrange,
                                 GtTimer *timer, GtError *err);

/* Returns true if <hcr_enc> was initialized with <descs> = true. */
bool          gt_hcr_encoder_has_descs_support(GtHcrEncoder *hcr_enc);

/* Sets the applied sampling method to <method> for <hcr_enc>. Default is
   pagewise. */
void          gt_hcr_encoder_set_sampling_none(GtHcrEncoder *hcr_enc);
void          gt_hcr_encoder_set_sampling_regular(GtHcrEncoder *hcr_enc);
void          gt_hcr_encoder_set_sampling_page(GtHcrEncoder *hcr_enc);

bool          gt_hcr_encoder_sampling_is_regular(GtHcrEncoder *hcr_enc);
bool          gt_hcr_encoder_sampling_is_page(GtHcrEncoder *hcr_enc);

/* Sets sampling rate to <srate> for the object <hcr_enc>. */
void          gt_hcr_encoder_set_sampling_rate(GtHcrEncoder *hcr_enc,
                                               unsigned long srate);

/* Returns the sampling rate of the object <hcr_enc>. */
unsigned long gt_hcr_encoder_get_sampling_rate(GtHcrEncoder *hcr_enc);

/* Encodes <hcr_enc> and writes the encoding to a file with base name <name>. */
int           gt_hcr_encoder_encode(GtHcrEncoder *hcr_enc, const char *name,
                                    GtTimer *timer, GtError *err);

/* Returns a new GtHcrDecoder object. <name> is the base name of a hcr encoded
   file. The description lines will not be decoded when <descs> is set to
   false.  */
GtHcrDecoder* gt_hcr_decoder_new(const char *name, GtAlphabet *alpha,
                                 bool descs, GtTimer *timer, GtError *err);

/* Returns true if <hcr_dec> was initialized with <descs> = true. */
bool          gt_hcr_decoder_has_descs_support(GtHcrDecoder *hcr_dec);

/* Decodes read with number <readnum> and writes decoding to the three char
   pointers. The addresses of the pointer must be allocated with sufficient
   memory space. <desc> gets reset and filled. */
int           gt_hcr_decoder_decode(GtHcrDecoder *hcr_dec,
                                    unsigned long readnum, char *seq,
                                    char *qual, GtStr * desc, GtError *err);

/* Decodes the hcr encoded file starting at record number <start> until record
   number <end> and writes the decoding to a file with base name <name>. */
int           gt_hcr_decoder_decode_range(GtHcrDecoder *hcr_dec,
                                          const char *name, unsigned long start,
                                          unsigned long end, GtTimer *timer,
                                          GtError *err);

/* Returns the total number of reads in <hcr_dec>. */
unsigned long gt_hcr_decoder_num_of_reads(GtHcrDecoder *hcr_dec);

/* Returns the read length of the reads in file with filenumber <filenum>. */
unsigned long gt_hcr_decoder_readlength(GtHcrDecoder *hcr_dec,
                                        unsigned long filenum);

void          gt_hcr_encoder_delete(GtHcrEncoder *hcr_enc);

void          gt_hcr_decoder_delete(GtHcrDecoder *hcr_dec);
#endif
