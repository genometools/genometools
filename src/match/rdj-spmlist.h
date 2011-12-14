/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_SPMLIST_H
#define RDJ_SPMLIST_H

#include "core/error_api.h"
#include "match/rdj-spmproc.h"

                 /* header byte: */
typedef enum {
  GT_SPMLIST_BIN32      = 2,
  GT_SPMLIST_BIN64      = 3,
  GT_SPMLIST_ASCII   /* = any other value */,
} GtSpmlistFormat;

#define DECLARE_GT_SPMLIST_BIN_FORMAT(BITS)\
void gt_spmlist_write_header_bin ## BITS(FILE *file);\
void gt_spmproc_show_bin ## BITS(unsigned long suffix_seqnum,\
    unsigned long prefix_seqnum, unsigned long length, bool suffixseq_direct,\
    bool prefixseq_direct, void *file)

DECLARE_GT_SPMLIST_BIN_FORMAT(32);
DECLARE_GT_SPMLIST_BIN_FORMAT(64);

/* parse a spmlist file; format is recognized by reading the first byte */
int gt_spmlist_parse(const char* filename, unsigned long min_length,
    GtSpmproc processoverlap, void *data, GtError *err);

void gt_spmproc_show_ascii(unsigned long suffix_seqnum,
    unsigned long prefix_seqnum, unsigned long length, bool suffixseq_direct,
    bool prefixseq_direct, void *data /* GtFile */);

/* approx overlaps: */

int gt_spmlist_approx_parse(const char* filename, unsigned long min_length,
    GtSpmprocA processoverlap, void *data, GtError *err);

void gt_spmproc_a_show_ascii(unsigned long suffix_seqnum,
    unsigned long prefix_seqnum, unsigned long suffix_length,
    unsigned long prefix_length, unsigned long unit_edist,
    bool suffixseq_direct, bool prefixseq_direct, void *data /* GtFile */);

/* tests */

int gt_spmlist_unit_test(GtError *err);

#endif
