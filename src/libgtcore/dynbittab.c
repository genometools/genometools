/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include <string.h>
#include "libgtcore/dynbittab.h"
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"

struct DynBittab {
  unsigned long *tabptr,
                tabsize,
                num_of_bits;
};

DynBittab* dynbittab_new(void)
{
  return ma_calloc(1, sizeof (DynBittab));
}

static unsigned long determine_tabsize(unsigned long num_of_bits)
{
  if (num_of_bits / (8UL * sizeof (unsigned long)))
    return 1 + ((num_of_bits - 1) / (8UL * sizeof (unsigned long)));
  return 1UL;
}

void dynbittab_set_bit(DynBittab *b, unsigned long bit)
{
  unsigned long new_tabsize;
  assert(b);
  /* make sure tab is large enough */
  if (bit >= b->num_of_bits) {
    if ((new_tabsize = determine_tabsize(bit + 1)) > b->tabsize) {
      b->tabptr = ma_realloc(b->tabptr, new_tabsize * sizeof (unsigned long));
      memset(b->tabptr + b->tabsize, 0,
             (new_tabsize - b->tabsize) * sizeof (unsigned long));
      b->tabsize = new_tabsize;
    }
    b->num_of_bits = bit + 1;
  }
  /* set bit */
  b->tabptr[(bit >> 3) / sizeof (unsigned long)] |=
    1UL << (bit & (8UL * sizeof (unsigned long) - 1));
}

void dynbittab_unset_bit(DynBittab *b, unsigned long bit)
{
  assert(b);
  if (bit < b->num_of_bits) {
    b->tabptr[(bit >> 3) / sizeof (unsigned long)] &=
      ~(1UL << (bit & (8UL * sizeof (unsigned long) - 1)));
  }
}

bool dynbittab_bit_is_set(const DynBittab *b, unsigned long bit)
{
  assert(b);
  if (bit < b->num_of_bits &&
      (b->tabptr[(bit >> 3) / sizeof (unsigned long)] &
       1UL << (bit & (8UL * sizeof (unsigned long) - 1)))) {
    return true;
  }
  return false;
}

int dynbittab_unit_test(Error *err)
{
  unsigned long i;
  DynBittab *b;
  int had_err = 0;

  b = dynbittab_new();
  for (i = 0; !had_err && i < 256; i++) {
    ensure(had_err, !dynbittab_bit_is_set(b, i));
  }
  if (!had_err) {
    dynbittab_set_bit(b, 0);
    dynbittab_set_bit(b, 32);
    dynbittab_set_bit(b, 64);
    dynbittab_set_bit(b, 77);
    dynbittab_set_bit(b, 96);
    dynbittab_set_bit(b, 123);
  }
  ensure(had_err, dynbittab_bit_is_set(b, 0));
  ensure(had_err, dynbittab_bit_is_set(b, 32));
  ensure(had_err, dynbittab_bit_is_set(b, 64));
  ensure(had_err, dynbittab_bit_is_set(b, 77));
  ensure(had_err, dynbittab_bit_is_set(b, 96));
  ensure(had_err, dynbittab_bit_is_set(b, 123));
  for (i = 124; !had_err && i < 256; i++) {
    ensure(had_err, !dynbittab_bit_is_set(b, i));
  }
  dynbittab_delete(b);

  b = dynbittab_new();
  for (i = 0; !had_err && i < 256; i++) {
    ensure(had_err, !dynbittab_bit_is_set(b, i));
  }
  if (!had_err) {
    dynbittab_set_bit(b, 1);
    dynbittab_set_bit(b, 33);
    dynbittab_set_bit(b, 65);
    dynbittab_set_bit(b, 77);
    dynbittab_set_bit(b, 97);
    dynbittab_set_bit(b, 124);
  }
  ensure(had_err, dynbittab_bit_is_set(b, 1));
  ensure(had_err, dynbittab_bit_is_set(b, 33));
  ensure(had_err, dynbittab_bit_is_set(b, 65));
  ensure(had_err, dynbittab_bit_is_set(b, 77));
  ensure(had_err, dynbittab_bit_is_set(b, 97));
  ensure(had_err, dynbittab_bit_is_set(b, 124));
  for (i = 125; !had_err && i < 256; i++) {
    ensure(had_err, !dynbittab_bit_is_set(b, i));
  }
  dynbittab_delete(b);

  return had_err;
}

void dynbittab_delete(DynBittab *b)
{
  if (!b) return;
  ma_free(b->tabptr);
  ma_free(b);
}
