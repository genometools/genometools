/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef QUALITY_H
#define QUALITY_H

#include "core/types_api.h"

/* Converts a FASTQ quality character in PHRED scale to PHRED quality score */
unsigned int gt_quality_fastq_to_phred(const GtUchar c);
/* Converts a FASTQ quality character in Solexa/Illumina scale to
   Solexa/Illumina quality score */
int          gt_quality_fastq_to_solexa(const GtUchar c);

/* Converts a PHRED quality score to Solexa/Illumina quality score */
double       gt_quality_phred_to_solexa(double phred);
/* Converts a Solexa/Illumina quality score to PHRED quality score */
double       gt_quality_solexa_to_phred(double solexascore);

/* Converts a PHRED quality score to probability that the base with that
   score was a calling error */
double       gt_quality_phred_to_errorprob(double phredscore);
/* Converts a Solexa/Illumina quality score to probability that the base with
   that score was a calling error */
double       gt_quality_solexa_to_errorprob(double solexascore);

int          gt_quality_unit_test(GtError *err);

#endif
