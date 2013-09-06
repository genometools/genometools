/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef STRING_MATCHING_H
#define STRING_MATCHING_H

#include <stdbool.h>

/* Gets called for every match. If it returns true, the matching is stopped. */
typedef bool (*GtProcessMatch)(GtUword pos, void *data);

/* Boyer-Moore-Horspool alg. (O(n*m) time worst-case; sublinear on average). */
void gt_string_matching_bmh(const char *s, GtUword n,
                            const char *p, GtUword m,
                            GtProcessMatch, void *data);

/* Knuth-Morris-Pratt algorithm (O(n+m) time).
   Returns the last common prefix length.*/
GtUword gt_string_matching_kmp(const char *s, GtUword n,
                                     const char *p, GtUword m,
                                     GtProcessMatch, void *data);

/* Shift-And algorithm (O(n*(m/|w|) time, |w| is the word size). */
void gt_string_matching_shift_and(const char *s, GtUword n,
                                  const char *p, GtUword m,
                                  GtProcessMatch, void *data);

/* Brute Force algorithm (O(n*m) time worst-case). */
void gt_string_matching_brute_force(const char *s, GtUword n,
                                    const char *p, GtUword m,
                                    GtProcessMatch, void *data);

int  gt_string_matching_unit_test(GtError*);

#endif
