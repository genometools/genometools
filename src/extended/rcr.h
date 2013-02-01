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

#ifndef RCR_H
#define RCR_H

#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/timer_api.h"

#define RCRFILESUFFIX ".rcr"

/* Classes <GtRcrEncoder> and <GtRcrDecoder> use mapped short reads stored as
   sam/bam and the corresponding reference sequences to compress these reads. */
typedef struct GtRcrEncoder GtRcrEncoder;
typedef struct GtRcrDecoder GtRcrDecoder;

/* Returns a new GtRcrEncoder object. <ref> points to a reference genome,
   <filename> is a BAM file containing alignments of short reads to
   the reference genome. If <vquals> is true, quality values of read
   positions having variations compared to the reference will be preserved.
   If <mquals> is true, the mapping quality of an alignment will be preserved.
   If <quals> is true, the quality values of all bases will be preserved.
   If <ureads> is true, unmapped reads will be written to a separated FASTQ.
   If <descs> is true, read names will be preserved. <vquals> and <quals>
   exclude each other. */
GtRcrEncoder* gt_rcr_encoder_new(const GtEncseq *ref,
                                 const char *filename,
                                 bool vquals,
                                 bool mquals,
                                 bool quals,
                                 bool ureads,
                                 bool descs,
                                 GtTimer *timer,
                                 GtError *err);

/* Enables verbosity for <rcr_enc>. That is, statistics about the encoding
   will be printed. */
void          gt_rcr_encoder_enable_verbosity(GtRcrEncoder *rcr_enc);

/* Disables verbosity for <rcr_enc>. */
void          gt_rcr_encoder_disable_verbosity(GtRcrEncoder *rcr_enc);

/* Writes the encoding of the BAM file associated with <rcr_enc> to a file
   given by <name> plus suffix ".rcr". */
int           gt_rcr_encoder_encode(GtRcrEncoder *rcr_enc,
                                    const char *name,
                                    GtTimer *timer,
                                    GtError *err);

/* Returns a new GtRcrDecoder object. <ref> points to a reference genome, <name>
   denotes the file path of a RCR file. */
GtRcrDecoder* gt_rcr_decoder_new(const char *name,
                                 const GtEncseq *ref,
                                 GtTimer *timer,
                                 GtError *err);

/* will load encoded descriptions and use these instead of just numbering the
   reads/alignments */
int           gt_rcr_decoder_enable_description_support(GtRcrDecoder *rcr_dec,
                                                        GtError *err);

/* Disables description support instead of names, the reads will get numbers */
void          gt_rcr_decoder_disable_description_support(GtRcrDecoder *rcr_dec);

/* Writes the decoding of the RCR file associated with <rcr_dec> to a file
   given by <name> plus suffix ".rcr.decoded" */
int           gt_rcr_decoder_decode(GtRcrDecoder *rcr_dec,
                                    const char *name,
                                    GtTimer *timer,
                                    GtError *err);

/* Deletes <rcr_enc>.*/
void          gt_rcr_encoder_delete(GtRcrEncoder *rcr_enc);

/* Deletes <rcr_dec>.*/
void          gt_rcr_decoder_delete(GtRcrDecoder *rcr_dec);

#endif
