/*
  Copyright (c) 2012 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef ENCDESC_H
#define ENCDESC_H

#include "core/error_api.h"
#include "core/str_array_api.h"
#include "core/timer_api.h"
#include "extended/cstr_iterator.h"

/* The <GtEncdesc> class stores a sequence description, e.g. a FASTA header,
   in a compressed form. This can save a lot of disk space or memory for
   repetitive headers, for example in multiple FASTA files with short reads. */
typedef struct GtEncdesc GtEncdesc;

/* The <GtEncdescEncoder> class can be used to encode FASTA header, delivering a
   <GtEncdesc>. */
typedef struct GtEncdescEncoder GtEncdescEncoder;

/* Returns a new <GtEncdescEncoder> object. */
GtEncdescEncoder* gt_encdesc_encoder_new(void);

/* Sets <timer> as the timer to be used in <ee> during the encoding process. */
void              gt_encdesc_encoder_set_timer(GtEncdescEncoder *ee,
                                               GtTimer *timer);

/* Returns the <GtTimer> set in <ee>, might be NULL. */
GtTimer*          gt_encdesc_encoder_get_timer(GtEncdescEncoder *ee);

/* Set the sampling method of <ee> to __none__. Sampling increases encoded size
   and decreases time for random access. */
void              gt_encdesc_encoder_set_sampling_none(GtEncdescEncoder *ee);
/* Set the sampling method of <ee> to either __page__wise sampling. Sampling
   increases encoded size and decreases time for random access. */
void              gt_encdesc_encoder_set_sampling_page(GtEncdescEncoder *ee);
/* Set the sampling method of <ee> to either __regular__ sampling. Sampling
   increases encoded size and decreases time for random access. */
void              gt_encdesc_encoder_set_sampling_regular(GtEncdescEncoder *ee);

/* Returns true if __page__wise sampling is set in <ee>. */
bool              gt_encdesc_encoder_sampling_is_page(GtEncdescEncoder *ee);
/* Returns true if __regular__ sampling is set in <ee>. */
bool              gt_encdesc_encoder_sampling_is_regular(GtEncdescEncoder *ee);

/* Sets the sampling rate */
void              gt_encdesc_encoder_set_sampling_rate(GtEncdescEncoder *ee,
                                                       GtUword sampling_rate);

/* Returns the samplingrate set in <ee>. Returns GT_UNDEF_UWORD if sampling is
   disabled. */
GtUword           gt_encdesc_encoder_get_sampling_rate(GtEncdescEncoder *ee);

/* Uses the settings in <ee> to encode the strings provided by <cstr_iterator>
   and writes them to a file with prefix <name>. Returns 0 on success, otherwise
   <err> is set accordingly. */
int               gt_encdesc_encoder_encode(GtEncdescEncoder *ee,
                                            GtCstrIterator *cstr_iterator,
                                            const char *name,
                                            GtError *err);

/* Loads a <GtEncdesc> from file with prefix <name> */
GtEncdesc*        gt_encdesc_load(const char *name, GtError *err);

/* Returns the number of encoded headers in <encdesc> */
GtUword           gt_encdesc_num_of_descriptions(GtEncdesc *encdesc);

/* Decodes description with number <num> and writes it to <desc>, which will be
   reset before writing to it. Returns 1 on success, 0 on EOF and -1 on
   error. <err> is set accordingly. */
int               gt_encdesc_decode(GtEncdesc *encdesc,
                                    GtUword num,
                                    GtStr *desc,
                                    GtError *err);

void              gt_encdesc_delete(GtEncdesc *encdesc);

void              gt_encdesc_encoder_delete(GtEncdescEncoder *ee);

int               gt_encdesc_unit_test(GtError *err);

#endif
