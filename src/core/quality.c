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

#include <math.h>
#include <ctype.h>
#include "core/assert_api.h"
#include "core/ensure.h"
#include "core/quality.h"
#include "core/mathsupport.h"

unsigned int gt_quality_fastq_to_phred(const GtUchar c)
{
  gt_assert(isprint((int) c));
  return ((unsigned int) c) - 33;
}

int gt_quality_fastq_to_solexa(const GtUchar c)
{
  gt_assert(isprint((int) c));
  return ((int) c) - 64;
}

double gt_quality_phred_to_solexa(double phredscore)
{
  return 10.0*log10(pow(10.0, ((double) phredscore/10.0)) - 1);
}

double gt_quality_solexa_to_phred(double solexascore)
{
  return 10.0*log10(pow(10.0, ((double) solexascore/10.0)) + 1);
}

double gt_quality_phred_to_errorprob(double phredscore)
{
  double rval;
  rval = pow(10.0, (double)(-0.1*phredscore));
  gt_assert(rval >= 0.0 && rval <= 1.0);
  return rval;
}

double gt_quality_solexa_to_errorprob(double solexascore)
{
  double rval;
  rval = 1.0/(1.0 + pow(10.0, (solexascore/10.0)));
  gt_assert(rval >= 0.0 && rval <= 1.0);
  return rval;
}

int gt_quality_unit_test(GtError *err)
{
  int had_err = 0;
  double val, phredscore, solexascore;

  gt_ensure(had_err, gt_quality_fastq_to_phred('!') == 0);
  gt_ensure(had_err, gt_quality_fastq_to_phred('{') == 90);
  val = gt_quality_phred_to_errorprob(gt_quality_fastq_to_phred('!'));
  gt_ensure(had_err, gt_double_equals_double(val, 1.0));
  val = gt_quality_phred_to_errorprob(gt_quality_fastq_to_phred('{'));
  gt_ensure(had_err, gt_double_equals_double(val, 0.000000001));

  gt_ensure(had_err, gt_quality_fastq_to_solexa('@') == 0);
  gt_ensure(had_err, gt_quality_fastq_to_solexa('!') == -31);
  gt_ensure(had_err, gt_quality_fastq_to_solexa('8') == -8);
  gt_ensure(had_err, gt_quality_fastq_to_solexa('{') == 59);

  phredscore = gt_quality_fastq_to_phred('A');
  solexascore = gt_quality_fastq_to_solexa('A'+31);
  gt_ensure(had_err, phredscore == solexascore);

  /* testing conversion (has rounding errors!)
     taken from Biopython test cases
     (http://portal.open-bio.org/pipermail/biopython-dev
             /2009-February/005385.html) */
  solexascore = gt_quality_phred_to_solexa(90);
  gt_ensure(had_err, gt_double_equals_double(solexascore, 89.999999995657035));
  solexascore = gt_quality_phred_to_solexa(50);
  gt_ensure(had_err, gt_double_equals_double(solexascore, 49.99995657033466));
  solexascore = gt_quality_phred_to_solexa(10);
  gt_ensure(had_err, gt_double_equals_double(solexascore, 9.5424250943932485));
  solexascore = gt_quality_phred_to_solexa(1);
  gt_ensure(had_err, gt_double_equals_double(solexascore, -5.8682532438011537));
  solexascore = gt_quality_phred_to_solexa(0.1);
  gt_ensure(had_err, gt_double_equals_double(solexascore, -16.32774717238372));

  phredscore = gt_quality_solexa_to_phred(90);
  gt_ensure(had_err, gt_double_equals_double(phredscore, 90.000000004342922));
  phredscore = gt_quality_solexa_to_phred(10);
  gt_ensure(had_err, gt_double_equals_double(phredscore, 10.41392685158225));
  phredscore = gt_quality_solexa_to_phred(0);
  gt_ensure(had_err, gt_double_equals_double(phredscore, 3.0102999566398116));
  phredscore = gt_quality_solexa_to_phred(-20);
  gt_ensure(had_err, gt_double_equals_double(phredscore, 0.043213737826425784));

  solexascore = gt_quality_fastq_to_solexa('!');
  phredscore = gt_quality_solexa_to_phred(solexascore);
  gt_ensure(had_err,
            gt_double_equals_double(phredscore, 0.0034483543102526788));
  solexascore = gt_quality_fastq_to_solexa('{');
  phredscore = gt_quality_solexa_to_phred(solexascore);
  gt_ensure(had_err, gt_double_equals_double(phredscore, 59.000005467440147));

  return had_err;
}
