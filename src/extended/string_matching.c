/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <limits.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/bittab.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/undef_api.h"
#include "extended/string_matching.h"

#define STRING_MATCHING_NUM_OF_TESTS        256
#define STRING_MATCHING_MAX_STRING_LENGTH   100000
#define STRING_MATCHING_MAX_PATTERN_LENGTH  66

void gt_string_matching_bmh(const char *s, unsigned long n,
                            const char *p, unsigned long m,
                            GtProcessMatch process_match, void *data)
{
  unsigned long i, j, pos, d[UCHAR_MAX];
  gt_assert(s && p);
  if (m > n || !m || !n) /* no match possible */
    return;
  /* preprocessing */
  for (i = 0; i < UCHAR_MAX; i++)
    d[i] = m;
  for (j = 0; j < m-1; j++)
    d[(unsigned) p[j]] = m-(j+1);
  /* searching */
  pos = 0;
  while (pos <= n-m) {
    j = m;
    while (j > 0 && s[pos+j-1] == p[j-1])
      j--;
    if (j == 0 && process_match) {
      if (process_match(pos, data))
        return;
    }
    pos += d[(unsigned) s[pos+m-1]];
  }
}

static unsigned long* compute_prefixtab(const char *p, unsigned long m)
{
  unsigned long i, vlen = 0, *prefixtab;
  char b;
  gt_assert(p);
  prefixtab = gt_malloc(sizeof (unsigned long) * (m+1));
  prefixtab[0] = GT_UNDEF_ULONG; /* paranoia */
  if (m)
    prefixtab[1] = 0;
  for (i = 2; i <= m; i++) {
    b = p[i-1];
    while (vlen > 0 && p[vlen] != b)
      vlen = prefixtab[vlen];
    if (p[vlen] == b)
      vlen++;
    prefixtab[i] = vlen;
  }
  return prefixtab;
}

unsigned long gt_string_matching_kmp(const char *s, unsigned long n,
                                     const char *p, unsigned long m,
                                     GtProcessMatch process_match, void *data)
{
  unsigned long *prefixtab,
                j = 0,   /* position in s corresponding to the first character
                            in p */
                cpl = 0; /* length of common prefix of s[j]..s[n-1] and p */
  char b, c;
  gt_assert(s && p);
  if (!m || !n) /* no match possible */
    return 0;
  prefixtab = compute_prefixtab(p, m);
  while (j + cpl < n) {
    if (cpl == m) {                     /* case (1)  */
      if (process_match) {
        if (process_match(j, data))
            break;
      }
      j = j + cpl - prefixtab[cpl];
      cpl = prefixtab[cpl];
    }
    else {                              /* case (2)  */
      b = s[j+cpl];
      c = p[cpl];
      if (b != c) {
        if (cpl > 0) {
          j = j + cpl - prefixtab[cpl]; /* case (2a) */
          cpl = prefixtab[cpl];
        }
        else
          j++;                          /* case (2b) */

      }
      else
        cpl++;                          /* case (2c) */
    }
  }
  /* do not miss match at the last possible position */
  if (cpl == m && process_match)
    process_match(j, data);
  gt_free(prefixtab);
  return cpl;
}

void gt_string_matching_shift_and(const char *s, unsigned long n,
                                  const char *p, unsigned long m,
                                  GtProcessMatch process_match, void *data)
{
  GtBittab *D, *B[UCHAR_MAX] = { NULL };
  unsigned long i, j;
  gt_assert(s && p);
  if (m > n || !m || !n) /* no match possible */
    return;
  /* preprocessing */
  for (j = 0; j < m; j++) {
    if (!B[(unsigned) p[j]])
      B[(unsigned) p[j]] = gt_bittab_new(m);
    gt_bittab_set_bit(B[(unsigned) p[j]], j);
  }
  /* searching */
  D = gt_bittab_new(m);
  for (i = 0; i < n; i++) {
    gt_bittab_shift_left_equal(D);
    gt_bittab_set_bit(D, 0);
    if (B[(unsigned) s[i]])
      gt_bittab_and_equal(D, B[(unsigned) s[i]]);
    else
      gt_bittab_unset(D);
    if (gt_bittab_bit_is_set(D, m - 1) && process_match) {
      if (process_match(i - m + 1, data))
        break;
    }
  }
  /* free */
  for (i = 0; i < UCHAR_MAX; i++)
    gt_bittab_delete(B[i]);
  gt_bittab_delete(D);
}

void gt_string_matching_brute_force(const char *s, unsigned long n,
                                    const char *p, unsigned long m,
                                    GtProcessMatch process_match, void *data)
{
  unsigned long i;
  gt_assert(s && p);
  if (m > n || !m || !n) /* no match possible */
    return;
  for (i = 0; i <= n - m; i++) {
    unsigned long j = 0;
    while (j < m && s[i+j] == p[j])
      j++;
    if (j == m && process_match) {
      if (process_match(i, data))
        break;
    }
  }
}

static bool store_first_match(unsigned long pos, void *data)
{
  unsigned long *match = data;
  gt_assert(match);
  *match = pos;
  return true;
}

static bool store_match(unsigned long pos, void *data)
{
  GtArray *positions = data;
  gt_assert(positions);
  gt_array_add(positions, pos);
  return false;
}

int gt_string_matching_unit_test(GtError *err)
{
  char s[STRING_MATCHING_MAX_STRING_LENGTH+1],
       p[STRING_MATCHING_MAX_PATTERN_LENGTH+1], *text = "foo";
  GtArray *brute_force_matches,
        *bmh_matches,
        *kmp_matches,
        *shift_and_matches;
  unsigned long i, brute_force_match, bmh_match, kmp_match, shift_and_match;
  int had_err = 0;

  gt_error_check(err);

  brute_force_matches = gt_array_new(sizeof (unsigned long));
  bmh_matches = gt_array_new(sizeof (unsigned long));
  kmp_matches = gt_array_new(sizeof (unsigned long));
  shift_and_matches = gt_array_new(sizeof (unsigned long));

  /* match the empty pattern */
  gt_string_matching_brute_force(text, strlen(text), "", 0, store_match,
                              brute_force_matches);
  gt_string_matching_bmh(text, strlen(text), "", 0, store_match, bmh_matches);
  gt_string_matching_kmp(text, strlen(text), "", 0, store_match, kmp_matches);
  gt_string_matching_shift_and(text, strlen(text), "", 0, store_match,
                            shift_and_matches);

  gt_ensure(had_err, !gt_array_size(brute_force_matches));
  gt_ensure(had_err, !gt_array_size(bmh_matches));
  gt_ensure(had_err, !gt_array_size(kmp_matches));
  gt_ensure(had_err, !gt_array_size(shift_and_matches));

  for (i = 0; !had_err && i < STRING_MATCHING_NUM_OF_TESTS; i++) {
    unsigned long j, n, m;
    /* generate random string and pattern */
    n = gt_rand_max(STRING_MATCHING_MAX_STRING_LENGTH);
    m = gt_rand_max(STRING_MATCHING_MAX_PATTERN_LENGTH);
    for (j = 0; j < n; j++)
      s[j] = gt_rand_char();
    s[n] = '\0';
    for (j = 0; j < m; j++)
      p[j] = gt_rand_char();
    p[m] = '\0';
    /* matching (first match) */
    brute_force_match = GT_UNDEF_ULONG;
    bmh_match = GT_UNDEF_ULONG;
    kmp_match = GT_UNDEF_ULONG;
    shift_and_match = GT_UNDEF_ULONG;
    gt_string_matching_brute_force(s, n, p, m, store_first_match,
                                &brute_force_match);
    gt_string_matching_bmh(s, n, p, m, store_first_match, &bmh_match);
    gt_string_matching_kmp(s, n, p, m, store_first_match, &kmp_match);
    gt_string_matching_shift_and(s, n, p, m, store_first_match,
                                 &shift_and_match);
    /* comparing (first match) */
    gt_ensure(had_err, brute_force_match == bmh_match);
    gt_ensure(had_err, brute_force_match == kmp_match);
    gt_ensure(had_err, brute_force_match == shift_and_match);
    /* matching (all matches) */
    gt_string_matching_brute_force(s, n, p, m, store_match,
                                   brute_force_matches);
    gt_string_matching_bmh(s, n, p, m, store_match, bmh_matches);
    gt_string_matching_kmp(s, n, p, m, store_match, kmp_matches);
    gt_string_matching_shift_and(s, n, p, m, store_match, shift_and_matches);
    /* comparing (all matches) */
    gt_ensure(had_err, gt_array_size(brute_force_matches) ==
                    gt_array_size(bmh_matches));
    gt_ensure(had_err, gt_array_size(brute_force_matches) ==
                    gt_array_size(kmp_matches));
    gt_ensure(had_err, gt_array_size(brute_force_matches) ==
                    gt_array_size(shift_and_matches));
    gt_ensure(had_err, !gt_array_cmp(brute_force_matches, bmh_matches));
    gt_ensure(had_err, !gt_array_cmp(brute_force_matches, kmp_matches));
    gt_ensure(had_err, !gt_array_cmp(brute_force_matches, shift_and_matches));
    /* reset */
    gt_array_reset(brute_force_matches);
    gt_array_reset(bmh_matches);
    gt_array_reset(kmp_matches);
    gt_array_reset(shift_and_matches);
  }

  gt_array_delete(shift_and_matches);
  gt_array_delete(bmh_matches);
  gt_array_delete(kmp_matches);
  gt_array_delete(brute_force_matches);

  return had_err;
}
