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
#include "libgtcore/bittab.h"
#include "libgtcore/ensure.h"
#include "libgtcore/fa.h"
#include "libgtcore/ma.h"
#include "libgtcore/mathsupport.h"
#include "libgtcore/undef.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"

#define NUM_OF_TESTS    50
#define MAX_SIZE        1024

struct Bittab {
  unsigned long *tabptr,
                tabsize,
                num_of_bits;
};

Bittab* bittab_new(unsigned long num_of_bits)
{
  Bittab *b;

  assert(num_of_bits);

  b = ma_malloc(sizeof (Bittab));
  b->num_of_bits = num_of_bits;

  if (num_of_bits / (8UL * sizeof (unsigned long)))
    b->tabsize = 1 + ((num_of_bits - 1) / (8UL * sizeof (unsigned long)));
  else
    b->tabsize = 1UL;

  b->tabptr = ma_calloc(b->tabsize, sizeof (unsigned long));

  return b;
}

void bittab_set_bit(Bittab *b, unsigned long bit)
{
  assert(b && bit < b->num_of_bits);
  b->tabptr[(bit >> 3) / sizeof (unsigned long)] |=
    1UL << (bit & (8UL * sizeof (unsigned long) - 1));
}

void bittab_unset_bit(Bittab *b, unsigned long bit)
{
  assert(b && bit < b->num_of_bits);
  b->tabptr[(bit >> 3) / sizeof (unsigned long)] &=
    ~(1UL << (bit & (8UL * sizeof (unsigned long) - 1)));
}

void bittab_complement(Bittab *dest, const Bittab *src)
{
  unsigned long i;

  assert(dest && src && dest->num_of_bits == src->num_of_bits);

  for (i = 0; i < src->tabsize - 1; i++)
    dest->tabptr[i] = ~src->tabptr[i];

  /* the ``last'' bittab gets special treatment to prevent that unused bits
     become set. this could disturb subsequent bittab_count_set_bits() calls. */
  dest->tabptr[src->tabsize - 1] = ~src->tabptr[src->tabsize - 1] &
                                   (~0UL >> (- src->num_of_bits +
                                             src->tabsize * 8UL *
                                             sizeof (unsigned long)));
}

void bittab_equal(Bittab *dest, const Bittab *src)
{
  unsigned long i;
  assert(dest && src && dest->num_of_bits == src->num_of_bits);
  for (i = 0; i < src->tabsize; i++)
    dest->tabptr[i] = src->tabptr[i];
}

void bittab_and(Bittab *dest, const Bittab *src1, const Bittab *src2)
{
  unsigned long i;
  assert(dest && src1 && src2);
  assert(dest->num_of_bits == src1->num_of_bits);
  assert(dest->num_of_bits == src2->num_of_bits);
  for (i = 0; i < src1->tabsize; i++)
    dest->tabptr[i] = src1->tabptr[i] & src2->tabptr[i];
}

void bittab_or(Bittab *dest, const Bittab *src1, const Bittab *src2)
{
  unsigned long i;
  assert(dest && src1 && src2);
  assert(dest->num_of_bits == src1->num_of_bits);
  assert(dest->num_of_bits == src2->num_of_bits);
  for (i = 0; i < src1->tabsize; i++)
    dest->tabptr[i] = src1->tabptr[i] | src2->tabptr[i];
}

void bittab_nand(Bittab *dest,
                  const Bittab *minuend,
                  const Bittab *subtrahend)
{
  unsigned long i;
  assert(dest && minuend && subtrahend);
  assert(dest->num_of_bits == minuend->num_of_bits);
  assert(minuend->num_of_bits == subtrahend->num_of_bits);
  for (i = 0; i < dest->tabsize; i++)
    dest->tabptr[i] = minuend->tabptr[i] & ~subtrahend->tabptr[i];
}

void bittab_and_equal(Bittab *dest, const Bittab *src)
{
  unsigned long i;
  assert(dest && src);
  assert(dest->num_of_bits == src->num_of_bits);
  for (i = 0; i < dest->tabsize; i++)
    dest->tabptr[i] &= src->tabptr[i];
}

void bittab_or_equal(Bittab *dest, const Bittab *src)
{
  unsigned long i;
  assert(dest && src);
  assert(dest->num_of_bits == src->num_of_bits);
  for (i = 0; i < dest->tabsize; i++)
    dest->tabptr[i] |= src->tabptr[i];
}

void bittab_shift_left_equal(Bittab *b)
{
  unsigned long i, new_carry, old_carry = 0;
  assert(b);
  for (i = 0; i < b->tabsize; i++) {
    new_carry = b->tabptr[i] & (1UL << (8UL * sizeof (unsigned long) - 1));
    b->tabptr[i] = (b->tabptr[i] << 1) | old_carry;
    old_carry = new_carry >> (8UL * sizeof (unsigned long) - 1);
  }
}

void bittab_shift_right_equal(Bittab *b)
{
  unsigned long i, new_carry, old_carry = 0;
  assert(b);
  for (i = b->tabsize; i > 0; i--) {
    new_carry = b->tabptr[i-1] & 1UL;
    b->tabptr[i-1] = (b->tabptr[i-1] >> 1) | old_carry;
    old_carry = new_carry << (8UL * sizeof (unsigned long) - 1);
  }
}

void bittab_unset(Bittab *b)
{
  unsigned long i;
  assert(b);
  for (i = 0; i < b->tabsize; i++)
    b->tabptr[i] = 0;
}

void bittab_get_all_bitnums(const Bittab *b, Array *bitnums)
{
  unsigned long i;
  assert(b && bitnums);
  for (i = 0; i < b->num_of_bits; i++)
    if (bittab_bit_is_set(b, i)) array_add(bitnums, i);
}

bool bittab_bit_is_set(const Bittab *b, unsigned long bit)
{
  assert(b && bit < b->num_of_bits);
  if (b->tabptr[(bit >> 3) / sizeof (unsigned long)] &
      1UL << (bit & (8UL * sizeof (unsigned long) - 1))) {
    return true;
  }
  return false;
}

bool bittab_is_true(const Bittab *b)
{
  unsigned long i;
  assert(b);
  for (i = 0; i < b->tabsize; i++) {
    if (b->tabptr[i])
      return true;
  }
  return false;
}

bool bittab_cmp(const Bittab *b1, const Bittab *b2)
{
  unsigned long i;
  assert(b1 && b2 && b1->num_of_bits == b2->num_of_bits);
  for (i = 0; i < b1->tabsize; i++) {
    if (b1->tabptr[i] != b2->tabptr[i])
      return false;
  }
  return true;
}

unsigned long bittab_size(Bittab *b)
{
  assert(b);
  return b->num_of_bits;
}

unsigned long bittab_get_first_bitnum(const Bittab *b)
{
  unsigned long i, rval = UNDEF_ULONG;
  assert(b);
  for (i = 0; i < b->num_of_bits; i++)
    if (bittab_bit_is_set(b, i)) {
      rval = i;
      break;
    }
  if (rval == UNDEF_ULONG)
    return b->num_of_bits;
  return rval;
}

unsigned long bittab_get_last_bitnum(const Bittab *b)
{
  assert(b);
  return b->num_of_bits;
}

unsigned long bittab_get_next_bitnum(const Bittab *b, unsigned long curnum)
{
  unsigned long i, rval = UNDEF_ULONG;

  assert(b);
  assert(curnum < b->num_of_bits);
  for (i = curnum + 1; i < b->num_of_bits; i++)
    if (bittab_bit_is_set(b, i)) {
      rval = i;
      break;
    }
  if (rval == UNDEF_ULONG)
    return b->num_of_bits;
  return rval;
}

unsigned long bittab_count_set_bits(Bittab *b)
{
  static const unsigned char bits_in_char[256] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2,
                                                   2, 3, 2, 3, 3, 4, 1, 2, 2, 3,
                                                   2, 3, 3, 4, 2, 3, 3, 4, 3, 4,
                                                   4, 5, 1, 2, 2, 3, 2, 3, 3, 4,
                                                   2, 3, 3, 4, 3, 4, 4, 5, 2, 3,
                                                   3, 4, 3, 4, 4, 5, 3, 4, 4, 5,
                                                   4, 5, 5, 6, 1, 2, 2, 3, 2, 3,
                                                   3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                                                   2, 3, 3, 4, 3, 4, 4, 5, 3, 4,
                                                   4, 5, 4, 5, 5, 6, 2, 3, 3, 4,
                                                   3, 4, 4, 5, 3, 4, 4, 5, 4, 5,
                                                   5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
                                                   4, 5, 5, 6, 5, 6, 6, 7, 1, 2,
                                                   2, 3, 2, 3, 3, 4, 2, 3, 3, 4,
                                                   3, 4, 4, 5, 2, 3, 3, 4, 3, 4,
                                                   4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                                   2, 3, 3, 4, 3, 4, 4, 5, 3, 4,
                                                   4, 5, 4, 5, 5, 6, 3, 4, 4, 5,
                                                   4, 5, 5, 6, 4, 5, 5, 6, 5, 6,
                                                   6, 7, 2, 3, 3, 4, 3, 4, 4, 5,
                                                   3, 4, 4, 5, 4, 5, 5, 6, 3, 4,
                                                   4, 5, 4, 5, 5, 6, 4, 5, 5, 6,
                                                   5, 6, 6, 7, 3, 4, 4, 5, 4, 5,
                                                   5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                                                   4, 5, 5, 6, 5, 6, 6, 7, 5, 6,
                                                   6, 7, 6, 7, 7, 8 };
  unsigned long i, j, counter = 0;
  assert(b);
  for (i = 0; i < b->tabsize; i++)
    for (j = 0; j < sizeof (unsigned long); j++)
      counter += bits_in_char[((b->tabptr[i] >> (j * 8)) & 0xffu)];
  return counter;
}

void bittab_show(const Bittab *b, FILE *outfp)
{
  unsigned long i;
  assert(b && outfp);
  /* header line */
  for (i = 0; i < b->num_of_bits; i++)
    fprintf(outfp, "%lu", i % 10);
  (void) putc('\n', outfp);
  /* actual bits */
  for (i = 0; i < b->num_of_bits; i++) {
    if (bittab_bit_is_set(b, i))
      (void) putc('1', outfp);
    else
      (void) putc('0', outfp);
  }
  (void) putc('\n', outfp);
}

int bittab_example(UNUSED Error *err)
{
  unsigned long bit;
  Bittab *b;
  error_check(err);

  b = bittab_new(32);
  bittab_set_bit(b, 8);
  bittab_set_bit(b, 16);
  bittab_set_bit(b, 24);

  /* a typical iterator loop */
  for (bit  = bittab_get_first_bitnum(b);
       bit != bittab_get_last_bitnum(b);
       bit  = bittab_get_next_bitnum(b, bit)) {
    /* ... */
  }

  bittab_delete(b);

  return 0;
}

int bittab_unit_test(Error *err)
{
  unsigned long i, j, size, bit, counter;
  Bittab *b, *tmp, *and;
  FILE *fp;
  int had_err = 0;
  error_check(err);

  for (i = 0; i < NUM_OF_TESTS && !had_err; i++) {
    size = rand_max(MAX_SIZE) + 1;
    b = bittab_new(size);
    tmp = bittab_new(size);
    and = bittab_new(size);
    ensure(had_err, bittab_size(b) == size);

    for (j = 0; j < size && !had_err; j++) {
      counter = 0;
      for (bit  = bittab_get_first_bitnum(b);
           bit != bittab_get_last_bitnum(b);
           bit  = bittab_get_next_bitnum(b, bit)) {
        counter++;
      }
      ensure(had_err, counter == j);

      ensure(had_err, bittab_count_set_bits(b) == j);
      ensure(had_err, !bittab_bit_is_set(b, j));
      bittab_set_bit(b, j);
      ensure(had_err, bittab_bit_is_set(b, j));

      bittab_complement(tmp, b);
      ensure(had_err, bittab_count_set_bits(tmp) == size - j - 1);
      ensure(had_err, !bittab_cmp(b, tmp));
      bittab_and(and, b, tmp);
      ensure(had_err, bittab_count_set_bits(and) == 0);

      bittab_unset(and);
      bittab_equal(and, b);
      bittab_or_equal(and, tmp);
      ensure(had_err, bittab_size(and) == size);

      bittab_equal(and, b);
      ensure(had_err, bittab_count_set_bits(and) == j + 1);
      bittab_and_equal(and, tmp);
      ensure(had_err, bittab_count_set_bits(and) == 0);

      bittab_complement(tmp, tmp);
      ensure(had_err, bittab_cmp(b, tmp));
    }

    ensure(had_err, bittab_count_set_bits(b) == size);
    bittab_complement(tmp, b);
    ensure(had_err, bittab_count_set_bits(tmp) == 0);

    for (j = 0; j < size && !had_err; j++) {
      bittab_unset_bit(b, j);
      ensure(had_err, !bittab_bit_is_set(b, j));
      ensure(had_err, bittab_count_set_bits(b) == size - j - 1);

      bittab_complement(tmp, b);
      ensure(had_err, !bittab_cmp(b, tmp));
      ensure(had_err, bittab_count_set_bits(tmp) == j + 1);
      bittab_and(and, b, tmp);
      ensure(had_err, bittab_count_set_bits(and) == 0);

      bittab_unset(and);
      bittab_equal(and, b);
      bittab_or_equal(and, tmp);
      ensure(had_err, bittab_size(and) == size);

      bittab_equal(and, b);
      ensure(had_err, bittab_count_set_bits(and) == size - j - 1);
      bittab_and_equal(and, tmp);
      ensure(had_err, bittab_count_set_bits(and) == 0);

      bittab_complement(tmp, tmp);
      ensure(had_err, bittab_cmp(b, tmp));
    }

    bittab_delete(b);
    bittab_delete(tmp);
    bittab_delete(and);
  }

  /* test bittab_show */
  fp = fa_xfopen("/dev/null", "w");
  b = bittab_new(80);
  for (i = 0; i < 80; i++) {
    if (i % 2)
      bittab_set_bit(b, i);
  }
  bittab_show(b, fp);
  bittab_delete(b);
  fa_xfclose(fp);

  /* test bittab_shift_left_equal() */
  b = bittab_new(125);
  bittab_set_bit(b, 0);
  bittab_set_bit(b, 32);
  bittab_set_bit(b, 64);
  bittab_set_bit(b, 77);
  bittab_set_bit(b, 96);
  bittab_set_bit(b, 123);
  ensure(had_err, bittab_count_set_bits(b) == 6);
  bittab_shift_left_equal(b);
  ensure(had_err, bittab_count_set_bits(b) == 6);
  ensure(had_err, bittab_bit_is_set(b, 1));
  ensure(had_err, bittab_bit_is_set(b, 33));
  ensure(had_err, bittab_bit_is_set(b, 65));
  ensure(had_err, bittab_bit_is_set(b, 78));
  ensure(had_err, bittab_bit_is_set(b, 97));
  ensure(had_err, bittab_bit_is_set(b, 124));
  bittab_delete(b);

  /* test bittab_shift_right_equal() */
  b = bittab_new(125);
  bittab_set_bit(b, 1);
  bittab_set_bit(b, 33);
  bittab_set_bit(b, 65);
  bittab_set_bit(b, 77);
  bittab_set_bit(b, 97);
  bittab_set_bit(b, 124);
  ensure(had_err, bittab_count_set_bits(b) == 6);
  bittab_shift_right_equal(b);
  ensure(had_err, bittab_count_set_bits(b) == 6);
  ensure(had_err, bittab_bit_is_set(b, 0));
  ensure(had_err, bittab_bit_is_set(b, 32));
  ensure(had_err, bittab_bit_is_set(b, 64));
  ensure(had_err, bittab_bit_is_set(b, 76));
  ensure(had_err, bittab_bit_is_set(b, 96));
  ensure(had_err, bittab_bit_is_set(b, 123));
  bittab_delete(b);

  return had_err;
}

void bittab_delete(Bittab *b)
{
  if (!b) return;
  ma_free(b->tabptr);
  ma_free(b);
}
