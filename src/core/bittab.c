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

#include "core/bittab.h"
#include "core/ensure.h"
#include "core/fa.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"

#define BITTAB_NUM_OF_TESTS  50
#define BITTAB_MAX_SIZE      1024

struct GtBittab {
  unsigned long *tabptr,
                tabsize,
                num_of_bits;
};

GtBittab* gt_bittab_new(unsigned long num_of_bits)
{
  GtBittab *b;

  gt_assert(num_of_bits != 0);

  b = gt_malloc(sizeof (GtBittab));
  b->num_of_bits = num_of_bits;

  if (num_of_bits / (8UL * sizeof (unsigned long)))
    b->tabsize = 1 + ((num_of_bits - 1) / (8UL * sizeof (unsigned long)));
  else
    b->tabsize = 1UL;

  b->tabptr = gt_calloc(b->tabsize, sizeof (unsigned long));

  return b;
}

void gt_bittab_set_bit(GtBittab *b, unsigned long bit)
{
  gt_assert(b && bit < b->num_of_bits);
  b->tabptr[(bit >> 3) / sizeof (unsigned long)] |=
    1UL << (bit & (8UL * sizeof (unsigned long) - 1));
}

void gt_bittab_unset_bit(GtBittab *b, unsigned long bit)
{
  gt_assert(b && bit < b->num_of_bits);
  b->tabptr[(bit >> 3) / sizeof (unsigned long)] &=
    ~(1UL << (bit & (8UL * sizeof (unsigned long) - 1)));
}

void gt_bittab_complement(GtBittab *dest, const GtBittab *src)
{
  unsigned long i;

  gt_assert(dest && src && dest->num_of_bits == src->num_of_bits);

  for (i = 0; i < src->tabsize - 1; i++)
    dest->tabptr[i] = ~src->tabptr[i];

  /* the ``last'' bittab gets special treatment to prevent that unused bits
     become set. this could disturb subsequent gt_bittab_count_set_bits() calls.
   */
  dest->tabptr[src->tabsize - 1] = ~src->tabptr[src->tabsize - 1] &
                                   (~0UL >> (- src->num_of_bits +
                                             src->tabsize * 8UL *
                                             sizeof (unsigned long)));
}

void gt_bittab_equal(GtBittab *dest, const GtBittab *src)
{
  unsigned long i;
  gt_assert(dest && src && dest->num_of_bits == src->num_of_bits);
  for (i = 0; i < src->tabsize; i++)
    dest->tabptr[i] = src->tabptr[i];
}

void gt_bittab_and(GtBittab *dest, const GtBittab *src1,
                   const GtBittab *src2)
{
  unsigned long i;
  gt_assert(dest && src1 && src2);
  gt_assert(dest->num_of_bits == src1->num_of_bits);
  gt_assert(dest->num_of_bits == src2->num_of_bits);
  for (i = 0; i < src1->tabsize; i++)
    dest->tabptr[i] = src1->tabptr[i] & src2->tabptr[i];
}

void gt_bittab_or(GtBittab *dest, const GtBittab *src1, const GtBittab *src2)
{
  unsigned long i;
  gt_assert(dest && src1 && src2);
  gt_assert(dest->num_of_bits == src1->num_of_bits);
  gt_assert(dest->num_of_bits == src2->num_of_bits);
  for (i = 0; i < src1->tabsize; i++)
    dest->tabptr[i] = src1->tabptr[i] | src2->tabptr[i];
}

void gt_bittab_nand(GtBittab *dest,
                  const GtBittab *minuend,
                  const GtBittab *subtrahend)
{
  unsigned long i;
  gt_assert(dest && minuend && subtrahend);
  gt_assert(dest->num_of_bits == minuend->num_of_bits);
  gt_assert(minuend->num_of_bits == subtrahend->num_of_bits);
  for (i = 0; i < dest->tabsize; i++)
    dest->tabptr[i] = minuend->tabptr[i] & ~subtrahend->tabptr[i];
}

void gt_bittab_and_equal(GtBittab *dest, const GtBittab *src)
{
  unsigned long i;
  gt_assert(dest && src);
  gt_assert(dest->num_of_bits == src->num_of_bits);
  for (i = 0; i < dest->tabsize; i++)
    dest->tabptr[i] &= src->tabptr[i];
}

void gt_bittab_or_equal(GtBittab *dest, const GtBittab *src)
{
  unsigned long i;
  gt_assert(dest && src);
  gt_assert(dest->num_of_bits == src->num_of_bits);
  for (i = 0; i < dest->tabsize; i++)
    dest->tabptr[i] |= src->tabptr[i];
}

void gt_bittab_shift_left_equal(GtBittab *b)
{
  unsigned long i, new_carry, old_carry = 0;
  gt_assert(b);
  for (i = 0; i < b->tabsize; i++) {
    new_carry = b->tabptr[i] & (1UL << (8UL * sizeof (unsigned long) - 1));
    b->tabptr[i] = (b->tabptr[i] << 1) | old_carry;
    old_carry = new_carry >> (8UL * sizeof (unsigned long) - 1);
  }
}

void gt_bittab_shift_right_equal(GtBittab *b)
{
  unsigned long i, new_carry, old_carry = 0;
  gt_assert(b);
  for (i = b->tabsize; i > 0; i--) {
    new_carry = b->tabptr[i-1] & 1UL;
    b->tabptr[i-1] = (b->tabptr[i-1] >> 1) | old_carry;
    old_carry = new_carry << (8UL * sizeof (unsigned long) - 1);
  }
}

void gt_bittab_unset(GtBittab *b)
{
  unsigned long i;
  gt_assert(b);
  for (i = 0; i < b->tabsize; i++)
    b->tabptr[i] = 0;
}

void gt_bittab_get_all_bitnums(const GtBittab *b, GtArray *bitnums)
{
  unsigned long i;
  gt_assert(b && bitnums);
  for (i = 0; i < b->num_of_bits; i++)
    if (gt_bittab_bit_is_set(b, i)) gt_array_add(bitnums, i);
}

bool gt_bittab_bit_is_set(const GtBittab *b, unsigned long bit)
{
  gt_assert(b && bit < b->num_of_bits);
  if (b->tabptr[(bit >> 3) / sizeof (unsigned long)] &
      1UL << (bit & (8UL * sizeof (unsigned long) - 1))) {
    return true;
  }
  return false;
}

bool gt_bittab_is_true(const GtBittab *b)
{
  unsigned long i;
  gt_assert(b);
  for (i = 0; i < b->tabsize; i++) {
    if (b->tabptr[i])
      return true;
  }
  return false;
}

bool gt_bittab_cmp(const GtBittab *b1, const GtBittab *b2)
{
  unsigned long i;
  gt_assert(b1 && b2 && b1->num_of_bits == b2->num_of_bits);
  for (i = 0; i < b1->tabsize; i++) {
    if (b1->tabptr[i] != b2->tabptr[i])
      return false;
  }
  return true;
}

unsigned long gt_bittab_size(GtBittab *b)
{
  gt_assert(b);
  return b->num_of_bits;
}

unsigned long gt_bittab_get_first_bitnum(const GtBittab *b)
{
  unsigned long i, rval = GT_UNDEF_ULONG;
  gt_assert(b);
  for (i = 0; i < b->num_of_bits; i++)
    if (gt_bittab_bit_is_set(b, i)) {
      rval = i;
      break;
    }
  if (rval == GT_UNDEF_ULONG)
    return b->num_of_bits;
  return rval;
}

unsigned long gt_bittab_get_last_bitnum(const GtBittab *b)
{
  gt_assert(b);
  return b->num_of_bits;
}

unsigned long gt_bittab_get_next_bitnum(const GtBittab *b,
                                        unsigned long curnum)
{
  unsigned long i, rval = GT_UNDEF_ULONG;

  gt_assert(b);
  gt_assert(curnum < b->num_of_bits);
  for (i = curnum + 1; i < b->num_of_bits; i++)
    if (gt_bittab_bit_is_set(b, i)) {
      rval = i;
      break;
    }
  if (rval == GT_UNDEF_ULONG)
    return b->num_of_bits;
  return rval;
}

unsigned long gt_bittab_count_set_bits(const GtBittab *b)
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
  gt_assert(b);
  for (i = 0; i < b->tabsize; i++)
    for (j = 0; j < sizeof (unsigned long); j++)
      counter += bits_in_char[((b->tabptr[i] >> (j * 8)) & 0xffu)];
  return counter;
}

void gt_bittab_show(const GtBittab *b, FILE *outfp)
{
  unsigned long i;
  gt_assert(b && outfp);
  /* header line */
  for (i = 0; i < b->num_of_bits; i++)
    fprintf(outfp, "%lu", i % 10);
  gt_xfputc('\n', outfp);
  /* actual bits */
  for (i = 0; i < b->num_of_bits; i++) {
    if (gt_bittab_bit_is_set(b, i))
      gt_xfputc('1', outfp);
    else
      gt_xfputc('0', outfp);
  }
  gt_xfputc('\n', outfp);
}

int gt_bittab_example(GT_UNUSED GtError *err)
{
  unsigned long bit;
  GtBittab *b;
  gt_error_check(err);

  b = gt_bittab_new(32);
  gt_bittab_set_bit(b, 8);
  gt_bittab_set_bit(b, 16);
  gt_bittab_set_bit(b, 24);

  /* a typical iterator loop */
  for (bit  = gt_bittab_get_first_bitnum(b);
       bit != gt_bittab_get_last_bitnum(b);
       bit  = gt_bittab_get_next_bitnum(b, bit)) {
    /* ... */
  }

  gt_bittab_delete(b);

  return 0;
}

int gt_bittab_unit_test(GtError *err)
{
  unsigned long i, j, size, bit, counter;
  GtBittab *b, *tmp, *and;
  FILE *fp;
  int had_err = 0;
  gt_error_check(err);

  for (i = 0; i < BITTAB_NUM_OF_TESTS && !had_err; i++) {
    size = gt_rand_max(BITTAB_MAX_SIZE) + 1;
    b = gt_bittab_new(size);
    tmp = gt_bittab_new(size);
    and = gt_bittab_new(size);
    gt_ensure(had_err, gt_bittab_size(b) == size);

    for (j = 0; j < size && !had_err; j++) {
      counter = 0;
      for (bit  = gt_bittab_get_first_bitnum(b);
           bit != gt_bittab_get_last_bitnum(b);
           bit  = gt_bittab_get_next_bitnum(b, bit)) {
        counter++;
      }
      gt_ensure(had_err, counter == j);

      gt_ensure(had_err, gt_bittab_count_set_bits(b) == j);
      gt_ensure(had_err, !gt_bittab_bit_is_set(b, j));
      gt_bittab_set_bit(b, j);
      gt_ensure(had_err, gt_bittab_bit_is_set(b, j));

      gt_bittab_complement(tmp, b);
      gt_ensure(had_err, gt_bittab_count_set_bits(tmp) == size - j - 1);
      gt_ensure(had_err, !gt_bittab_cmp(b, tmp));
      gt_bittab_and(and, b, tmp);
      gt_ensure(had_err, gt_bittab_count_set_bits(and) == 0);

      gt_bittab_unset(and);
      gt_bittab_equal(and, b);
      gt_bittab_or_equal(and, tmp);
      gt_ensure(had_err, gt_bittab_size(and) == size);

      gt_bittab_equal(and, b);
      gt_ensure(had_err, gt_bittab_count_set_bits(and) == j + 1);
      gt_bittab_and_equal(and, tmp);
      gt_ensure(had_err, gt_bittab_count_set_bits(and) == 0);

      gt_bittab_complement(tmp, tmp);
      gt_ensure(had_err, gt_bittab_cmp(b, tmp));
    }

    gt_ensure(had_err, gt_bittab_count_set_bits(b) == size);
    gt_bittab_complement(tmp, b);
    gt_ensure(had_err, gt_bittab_count_set_bits(tmp) == 0);

    for (j = 0; j < size && !had_err; j++) {
      gt_bittab_unset_bit(b, j);
      gt_ensure(had_err, !gt_bittab_bit_is_set(b, j));
      gt_ensure(had_err, gt_bittab_count_set_bits(b) == size - j - 1);

      gt_bittab_complement(tmp, b);
      gt_ensure(had_err, !gt_bittab_cmp(b, tmp));
      gt_ensure(had_err, gt_bittab_count_set_bits(tmp) == j + 1);
      gt_bittab_and(and, b, tmp);
      gt_ensure(had_err, gt_bittab_count_set_bits(and) == 0);

      gt_bittab_unset(and);
      gt_bittab_equal(and, b);
      gt_bittab_or_equal(and, tmp);
      gt_ensure(had_err, gt_bittab_size(and) == size);

      gt_bittab_equal(and, b);
      gt_ensure(had_err, gt_bittab_count_set_bits(and) == size - j - 1);
      gt_bittab_and_equal(and, tmp);
      gt_ensure(had_err, gt_bittab_count_set_bits(and) == 0);

      gt_bittab_complement(tmp, tmp);
      gt_ensure(had_err, gt_bittab_cmp(b, tmp));
    }

    gt_bittab_delete(b);
    gt_bittab_delete(tmp);
    gt_bittab_delete(and);
  }

  /* test gt_bittab_show */
  fp = gt_fa_xfopen("/dev/null", "w");
  b = gt_bittab_new(80);
  for (i = 0; i < 80; i++) {
    if (i % 2)
      gt_bittab_set_bit(b, i);
  }
  gt_bittab_show(b, fp);
  gt_bittab_delete(b);
  gt_fa_xfclose(fp);

  /* test gt_bittab_shift_left_equal() */
  b = gt_bittab_new(125);
  gt_bittab_set_bit(b, 0);
  gt_bittab_set_bit(b, 32);
  gt_bittab_set_bit(b, 64);
  gt_bittab_set_bit(b, 77);
  gt_bittab_set_bit(b, 96);
  gt_bittab_set_bit(b, 123);
  gt_ensure(had_err, gt_bittab_count_set_bits(b) == 6);
  gt_bittab_shift_left_equal(b);
  gt_ensure(had_err, gt_bittab_count_set_bits(b) == 6);
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 1));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 33));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 65));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 78));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 97));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 124));
  gt_bittab_delete(b);

  /* test gt_bittab_shift_right_equal() */
  b = gt_bittab_new(125);
  gt_bittab_set_bit(b, 1);
  gt_bittab_set_bit(b, 33);
  gt_bittab_set_bit(b, 65);
  gt_bittab_set_bit(b, 77);
  gt_bittab_set_bit(b, 97);
  gt_bittab_set_bit(b, 124);
  gt_ensure(had_err, gt_bittab_count_set_bits(b) == 6);
  gt_bittab_shift_right_equal(b);
  gt_ensure(had_err, gt_bittab_count_set_bits(b) == 6);
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 0));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 32));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 64));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 76));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 96));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 123));
  gt_bittab_delete(b);

  return had_err;
}

void gt_bittab_delete(GtBittab *b)
{
  if (!b) return;
  gt_free(b->tabptr);
  gt_free(b);
}
