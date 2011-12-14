/*
  Copyright (c) 2009-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2009-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_SPMPROC_H
#define RDJ_SPMPROC_H

#include <stdbool.h>
#include "core/intbits.h"            /* GtBitsequence */
#include "core/error_api.h"          /* GtError */

/* Prototype for functions processing exact/approximate overlaps */

/* Exact: */
typedef void(*GtSpmproc)(unsigned long /*suffix_seqnum*/,
  unsigned long /*prefix_seqnum*/, unsigned long /*length*/,
  bool /*suffixseq_direct*/, bool /*prefixseq_direct*/, void* /*data*/);

/* Approximate: */
typedef void(*GtSpmprocA)(unsigned long /*suffix_seqnum*/,
  unsigned long /*prefix_seqnum*/, unsigned long /*suffix_length*/,
  unsigned long /*prefix_length*/, unsigned long /*unit_edist*/,
  bool /*suffixseq_direct*/, bool /*prefixseq_direct*/, void* /*data*/);

/* structs containing an spmproc and its data */
typedef struct {GtSpmproc proc; void *data;} GtSpmprocWithData;
typedef struct {GtSpmprocA proc; void *data;} GtSpmprocAWithData;
typedef union {GtSpmprocWithData e; GtSpmprocAWithData a;} GtSpmprocXWithData;

/* allow to use GtSpmproc where GtSpmprocA is expected skipping overlaps with
   unit edist > 0; data must be a pointer to GtSpmprocWithData */
void gt_spmproc_a_e(unsigned long suffix_seqnum, unsigned long prefix_seqnum,
  unsigned long suffix_length, unsigned long prefix_length,
  unsigned long unit_edist, bool suffixseq_direct, bool prefixseq_direct,
  void *data);

void gt_spmproc_count(unsigned long suffix_seqnum, unsigned long prefix_seqnum,
  unsigned long length, bool suffixseq_direct, bool prefixseq_direct,
  void *data);

/*
Given a GtBitsequence (large as the number of reads), acts as a filter
and calls the Spmproc outproc if the bit for both seqnums is not set.
The void* data must be of type GtSpmprocSkipData.
*/

typedef struct {
  GtBitsequence *to_skip;
  GtSpmprocXWithData out;
  unsigned long skipped_counter;
} GtSpmprocSkipData;

void gt_spmproc_skip(unsigned long suffix_seqnum, unsigned long prefix_seqnum,
  unsigned long length, bool suffixseq_direct, bool prefixseq_direct,
  void *data);

void gt_spmproc_a_skip(unsigned long suffix_seqnum,
  unsigned long prefix_seqnum, unsigned long suffix_length,
  unsigned long prefix_length, unsigned long unit_edist, bool suffixseq_direct,
  bool prefixseq_direct, void *data);

int gt_spmproc_skip_unit_test(GtError *err);

#endif
