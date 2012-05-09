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

#ifndef RDJ_CONTFINDER_H
#define RDJ_CONTFINDER_H

#include <stdint.h>
#include "core/file.h"

typedef struct GtContfinder GtContfinder;

/* totallength: if 0, the sum of the size of the files is used */
GtContfinder* gt_contfinder_new(GtStrArray *filenames, GtStr *indexname,
    bool output_encseq, GtError *err);

typedef enum {
  GT_CONTFINDER_SEQNUMS,
  GT_CONTFINDER_2BIT,
  GT_CONTFINDER_FASTA,
  GT_CONTFINDER_QUIET
} GtContfinderOutputFormat;

/*
 * rev: true = reverse complements are also considered
 *
 * format: QUIET = no output; otherwise output to GtFile *outfp
 *         in the specified format (seqnums = newline separated list of seqnums)
 *
 * sorted: false = output non contained reads, in their input order
 *         true = output all reads (including contained), lexicogr. sorted
 *
 * cntlistfilename: NULL = do not output contained reads list
 *                  otherwise = output binary contained reads list to file
 * sepposfilename: NULL = do not output separator positions list
 *                 otherwise = output binary contained reads list to file
 *                 [note: the totallength is output as last value]
 */
int gt_contfinder_run(GtContfinder *contfinder, bool rev, GtFile *outfp,
    GtContfinderOutputFormat format, bool sorted, const char *cntlistfilename,
    const char *sepposfilename, const char *copynumfilename,
    bool output_encseq, GtError *err);

void gt_contfinder_delete(GtContfinder *contfinder);

unsigned long gt_contfinder_totallength_without_sep(GtContfinder *contfinder);

unsigned long gt_contfinder_nofseqs(GtContfinder *contfinder);

unsigned long gt_contfinder_nofcontained(GtContfinder *contfinder);

unsigned long gt_contfinder_nofdiscarded(GtContfinder *contfinder);
unsigned long gt_contfinder_discarded_length(GtContfinder *contfinder);

unsigned long gt_contfinder_read_length(GtContfinder *contfinder);

void gt_contfinder_radixsort_str_eqlen_tester(GtContfinder *contfinder,
    bool mirrored, unsigned long depth,
    unsigned long maxdepth, bool print);

#endif
