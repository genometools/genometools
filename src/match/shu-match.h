/*
  Copyright (c) 2010 Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
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

#ifndef SHU_MATCH_H
#define SHU_MATCH_H

#include "match/eis-bwtseq.h"

/* threshold used for calculation of divergence */
#define THRESHOLD pow(10, -9)
#define DEFAULT_E 1e-3
#define DEFAULT_T 1e-5
#define DEFAULT_M DBL_MIN
#define N_CHOOSE 1000
#define N_S 1000    /* at the moment - number of elements in the s1 array */

/* returns the length of the matching prefix +1, that is it returns the
 * length of the shortest absent prefix */
unsigned int gt_pck_getShuStringLength(const BWTSeq* bwtSubject,
                                       const GtUchar* query,
                                       size_t queryLength);

/* returns the ½gc-content of the subject */
double gt_pck_getGCcontent(const BWTSeq *bwtSubject,
                           const GtAlphabet *alphabet);

/* calculates the divergence from the shulength */
double divergence (double E, /* relative error for shulen length */
                   double T, /* absolute error */
                   double M,
                   double shulen, /*avg shulength*/
                   unsigned long seqLen, /*querylength*/
                   double gc, /*combined gc-content*/
                   int n_choose,
                   double **lnChoose);

/* Jukes-Cantor transform of divergence: Kr */
double calculateKr(double d);

/* gives a 2-dim array used for precalculations, just for speed */
double **initializeLnChoose(unsigned long dim);

/* frees the precalc-matrix */
void freeLnChoose(double **matrix,
                  unsigned long dim);
#endif
