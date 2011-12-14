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

#ifndef RDJ_OVLFIND_DP_H
#define RDJ_OVLFIND_DP_H

#include "core/error_api.h"               /* GtError       */
#include "match/rdj-contfind-def.h"       /* GtContfind    */
#include "match/rdj-pairwise.h"           /* GtOvlfindMode */

/* Dynamic programming based overlap finder */

/*
  known issues:
  - case sensitive (a != A)
  - no notion of wildcards (i.e. n == n)
*/

GtContfind gt_ovlfind_dp(const char *u, unsigned long m,
    const char *v /* use NULL for self-comparison of u */, unsigned long n,
    double max_error, GtOvlfindMode mode, unsigned long min_length,
    bool find_submaximal, void (*smpproc) (unsigned long /* length on u */,
    unsigned long /* length on v */, unsigned long /* unit edit distance */,
    bool /* true if suffix comes from u, false if suffix comes from v */,
    void* /* procdata */), void* smpprocdata);

int gt_ovlfind_dp_unit_test(GtError *err);

#endif
