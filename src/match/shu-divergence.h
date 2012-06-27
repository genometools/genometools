/*
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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

#ifndef SHU_DIVERGENCE_H
#define SHU_DIVERGENCE_H

/* calculates the divergence from the shulength */
double gt_divergence (double E, /* relative error for shulen length */
                   double T, /* absolute error */
                   double M,
                   double threshold, /* break criteria for pmax calculation */
                   double shulen, /*avg shulength*/
                   unsigned long subjectLength, /*subjectlength*/
                   double gc, /*combined gc-content*/
                   double *ln_n_fac,
                   unsigned long n_s /*length of s-array*/);

/* Jukes-Cantor transform of divergence: Kr */
double gt_calculateKr(double d);

/* Calculates ln(n!) for numbers 0 to given n */
double *gt_get_ln_n_fac(unsigned long n);

#endif
