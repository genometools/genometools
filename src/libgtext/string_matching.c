/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <assert.h>
#include <limits.h>
#include <string.h>
#include "libgtcore/bittab.h"
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/mathsupport.h"
#include "libgtcore/undef.h"
#include "libgtext/string_matching.h"

#define NUM_OF_TESTS        256
#define MAX_STRING_LENGTH   100000
#define MAX_PATTERN_LENGTH  66

void string_matching_bmh(const char *s, unsigned long n,
                         const char *p, unsigned long m,
                         ProcessMatch process_match, void *data)
{
  unsigned long i, j, pos, d[UCHAR_MAX];
  assert(s && p);
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
  assert(p);
  prefixtab = ma_malloc(sizeof (unsigned long) * (m+1));
  prefixtab[0] = UNDEF_ULONG; /* paranoia */
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

unsigned long string_matching_kmp(const char *s, unsigned long n,
                                  const char *p, unsigned long m,
                                  ProcessMatch process_match, void *data)
{
  unsigned long *prefixtab,
                j = 0,   /* position in s corresponding to the first character
                            in p */
                cpl = 0; /* length of common prefix of s[j]..s[n-1] and p */
  char b, c;
  assert(s && p);
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
  ma_free(prefixtab);
  return cpl;
}

void string_matching_shift_and(const char *s, unsigned long n,
                               const char *p, unsigned long m,
                               ProcessMatch process_match, void *data)
{
  Bittab *D, *B[UCHAR_MAX] = { NULL };
  unsigned long i, j;
  assert(s && p);
  if (m > n || !m || !n) /* no match possible */
    return;
  /* preprocessing */
  for (j = 0; j < m; j++) {
    if (!B[(unsigned) p[j]])
      B[(unsigned) p[j]] = bittab_new(m);
    bittab_set_bit(B[(unsigned) p[j]], j);
  }
  /* searching */
  D = bittab_new(m);
  for (i = 0; i < n; i++) {
    bittab_shift_left_equal(D);
    bittab_set_bit(D, 0);
    if (B[(unsigned) s[i]])
      bittab_and_equal(D, B[(unsigned) s[i]]);
    else
      bittab_unset(D);
    if (bittab_bit_is_set(D, m - 1) && process_match) {
      if (process_match(i - m + 1, data))
        break;
    }
  }
  /* free */
  for (i = 0; i < UCHAR_MAX; i++)
    bittab_delete(B[i]);
  bittab_delete(D);
}

void string_matching_brute_force(const char *s, unsigned long n,
                                 const char *p, unsigned long m,
                                 ProcessMatch process_match, void *data)
{
  unsigned long i;
  assert(s && p);
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
  assert(match);
  *match = pos;
  return true;
}

static bool store_match(unsigned long pos, void *data)
{
  Array *positions = data;
  assert(positions);
  array_add(positions, pos);
  return false;
}

int string_matching_unit_test(Error *err)
{
  char s[MAX_STRING_LENGTH+1], p[MAX_PATTERN_LENGTH+1], *text = "foo";
  Array *brute_force_matches,
        *bmh_matches,
        *kmp_matches,
        *shift_and_matches;
  unsigned long i, brute_force_match, bmh_match, kmp_match, shift_and_match;
  int had_err = 0;

  error_check(err);

  brute_force_matches = array_new(sizeof (unsigned long));
  bmh_matches = array_new(sizeof (unsigned long));
  kmp_matches = array_new(sizeof (unsigned long));
  shift_and_matches = array_new(sizeof (unsigned long));

  /* match the empty pattern */
  string_matching_brute_force(text, strlen(text), "", 0, store_match,
                              brute_force_matches);
  string_matching_bmh(text, strlen(text), "", 0, store_match, bmh_matches);
  string_matching_kmp(text, strlen(text), "", 0, store_match, kmp_matches);
  string_matching_shift_and(text, strlen(text), "", 0, store_match,
                            shift_and_matches);

  ensure(had_err, !array_size(brute_force_matches));
  ensure(had_err, !array_size(bmh_matches));
  ensure(had_err, !array_size(kmp_matches));
  ensure(had_err, !array_size(shift_and_matches));

  for (i = 0; !had_err && i < NUM_OF_TESTS; i++) {
    unsigned long j, n, m;
    /* generate random string and pattern */
    n = rand_max(MAX_STRING_LENGTH);
    m = rand_max(MAX_PATTERN_LENGTH);
    for (j = 0; j < n; j++)
      s[j] = rand_char();
    s[n] = '\0';
    for (j = 0; j < m; j++)
      p[j] = rand_char();
    p[m] = '\0';
    /* matching (first match) */
    brute_force_match = UNDEF_ULONG;
    bmh_match = UNDEF_ULONG;
    kmp_match = UNDEF_ULONG;
    shift_and_match = UNDEF_ULONG;
    string_matching_brute_force(s, n, p, m, store_first_match,
                                &brute_force_match);
    string_matching_bmh(s, n, p, m, store_first_match, &bmh_match);
    string_matching_kmp(s, n, p, m, store_first_match, &kmp_match);
    string_matching_shift_and(s, n, p, m, store_first_match, &shift_and_match);
    /* comparing (first match) */
    ensure(had_err, brute_force_match == bmh_match);
    ensure(had_err, brute_force_match == kmp_match);
    ensure(had_err, brute_force_match == shift_and_match);
    /* matching (all matches) */
    string_matching_brute_force(s, n, p, m, store_match, brute_force_matches);
    string_matching_bmh(s, n, p, m, store_match, bmh_matches);
    string_matching_kmp(s, n, p, m, store_match, kmp_matches);
    string_matching_shift_and(s, n, p, m, store_match, shift_and_matches);
    /* comparing (all matches) */
    ensure(had_err, array_size(brute_force_matches) == array_size(bmh_matches));
    ensure(had_err, array_size(brute_force_matches) == array_size(kmp_matches));
    ensure(had_err, array_size(brute_force_matches) ==
                    array_size(shift_and_matches));
    ensure(had_err, !array_cmp(brute_force_matches, bmh_matches));
    ensure(had_err, !array_cmp(brute_force_matches, kmp_matches));
    ensure(had_err, !array_cmp(brute_force_matches, shift_and_matches));
    /* reset */
    array_reset(brute_force_matches);
    array_reset(bmh_matches);
    array_reset(kmp_matches);
    array_reset(shift_and_matches);
  }

  array_delete(shift_and_matches);
  array_delete(bmh_matches);
  array_delete(kmp_matches);
  array_delete(brute_force_matches);

  return had_err;
}
