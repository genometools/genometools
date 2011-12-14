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

#ifndef RDJ_OVLFIND_KMP_H
#define RDJ_OVLFIND_KMP_H

#include "core/error_api.h"          /* GtError       */
#include "match/rdj-contfind-def.h"  /* GtContfind    */
#include "match/rdj-pairwise.h"      /* GtOvlfindMode */

/* KMP overlap finder */

/*
  - case sensitive (a != A)
  - no notion of wildcards (i.e. n == n)
*/

#include <stdint.h>
#include "core/error.h"

typedef uint16_t gt_kmp_t;
#define GT_KMP_MAX UINT16_MAX

gt_kmp_t* gt_kmp_preproc(const char *seq, unsigned long seqlen);

GtContfind gt_ovlfind_kmp(const char *u, unsigned long u_length,
                          const gt_kmp_t *u_pi,
                          const char *v /* use NULL for u vs. u */,
                          unsigned long v_length, const gt_kmp_t *v_pi,
                          GtOvlfindMode m, unsigned long min_length,
                          bool find_nonmaximal,
                          void(*spmproc)
                            (unsigned long /* overlap length */,
                             bool /* true if suffix of u == prefix of v
                                     false if prefix of u == suffix of v */,
                             void* /* spmprocdata */),
                          void* spmprocdata);

int gt_kmp_preproc_unit_test(GtError *err);
int gt_ovlfind_kmp_unit_test(GtError *err);

#endif
