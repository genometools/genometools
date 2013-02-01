/*
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

#ifndef SAMPLING_H
#define SAMPLING_H

#include <stdbool.h>
#include <stdio.h>
#ifndef S_SPLINT_S
#include <sys/types.h>
#endif

#define GT_SAMPLING_DEFAULT_REGULAR_RATE 10000UL
#define GT_SAMPLING_DEFAULT_PAGE_RATE 100UL

/* Class <GtSampling> can be used to collect sampling information for variable
   length data that is written to a file. For each sampled element, it stores
   the number of that element and its page offset in the file. When read from
   a file, it helps to find the page offset of a sample. */
typedef struct GtSampling GtSampling;

/* Returns a new <GtSampling> object which uses regular sampling, <rate> sets
   the sampling rate and <first_offset> is the position of the first sample in
   the file, has to be a multiple of the pagesize. */
GtSampling*   gt_sampling_new_regular(unsigned long rate, off_t first_offset);

/* Returns a new <GtSampling> object which uses page oriented sampling, that is:
   the sampling rate defines how many pages should be filled before sampling the
   next element. <rate> sets the sampling rate and <first_offset> is the
   position of the first sample in the file, has to be a multiple of the
   pagesize. */
GtSampling*   gt_sampling_new_page(unsigned long rate, off_t first_offset);

/* Writes <sampling> to <FILE> <fp>. */
void          gt_sampling_write(GtSampling *sampling, FILE *fp);

/* Reads from <fp> the information previously stored in <fp> by
   <gt_sampling_write> and returns a new <GtSampling> object. <fp> has to point
   to the correct start position in the file. */
GtSampling*   gt_sampling_read(FILE *fp);

/* sets <*sampled_element> to the largest sampled element <= <element_num>, and
   sets <*position> to the offset where that sample starts. Returns -1 in case
   of error. */
void          gt_sampling_get_page(GtSampling *sampling,
                                   unsigned long element_num,
                                   unsigned long *sampled_element,
                                   size_t *position);

/* Returns the number of the sampled element of the current sample, this changes
   when <gt_sampling_get_page> or <gt_sampling_get_next_sample> are called */
unsigned long gt_sampling_get_current_elementnum(GtSampling *sampling);

/* This does NOT iterate, it just reports the number of the next sampled
   element, returns 0 if current sample is the last sampled. */
unsigned long gt_sampling_get_next_elementnum(GtSampling *sampling);

/* Sets <*sampled_element> to the next start of sampled elements, and
   <*position> to the corresponding offset in the file. This is an iterator, so
   subsequent calls will have different results. Returns 1 if sample was found,
   0 if there is no next sample, -1 on error. */
int           gt_sampling_get_next_sample(GtSampling *sampling,
                                          unsigned long *sampled_element,
                                          size_t *position);

/* Returns the sampling rate of <sampling>. */
unsigned long gt_sampling_get_rate(GtSampling *sampling);

/* Returns true if <sampling> uses regular sampling. */
bool          gt_sampling_is_regular(GtSampling *sampling);

/* Tells <sampling> to store a new sample with element number <element_num>
   and file offset <position>, where position is expected to be a multiple of
   pagesize. */
void          gt_sampling_add_sample(GtSampling *sampling,
                                     size_t position,
                                     unsigned long element_num);

/* TODO: document me */
bool          gt_sampling_is_next_element_sample(
                                          GtSampling *sampling,
                                          unsigned long pages_written,
                                          unsigned long elements_written,
                                          unsigned long elem_bit_size,
                                          unsigned long free_pagespace_bitsize);

void          gt_sampling_delete(GtSampling *sampling);

#endif
