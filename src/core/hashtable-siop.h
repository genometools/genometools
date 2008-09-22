/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#ifndef HASHTABLE_SIOP_H
#define HASHTABLE_SIOP_H

#include <limits.h>
#include <string.h>

#include "core/unused_api.h"

#define GOLDEN_RATIO_MULTIPLIER 2654435761UL

static inline uint32_t
ht_rotate_left_u32(uint32_t m, unsigned short k)
{
  return (((m)<<(k)) | ((m)>>((sizeof (m) * CHAR_BIT)-(k))));
}

static inline uint64_t
ht_rotate_left_u64(uint64_t m, unsigned short k)
{
  return (((m)<<(k)) | ((m)>>((sizeof (m) * CHAR_BIT)-(k))));
}

static inline uint32_t
ht_rotate_right_u32(uint32_t m, unsigned short k)
{
  return (((m)>>(k)) | ((m)<<((sizeof (m) * CHAR_BIT)-(k))));
}

static inline uint64_t
ht_rotate_right_u64(uint64_t m, unsigned short k)
{
  return (((m)>>(k)) | ((m)<<((sizeof (m) * CHAR_BIT)-(k))));
}

/*
 * The finalize function and the multiplicative hash functions is
 * derived from Bob Jenkins' public domain implementation
 * (http://www.burtleburtle.net/bob/)
 */

static inline uint32_t
ht_finalize3_u32(uint32_t a, uint32_t b, uint32_t c)
{
  c ^= b; c -= ht_rotate_left_u32(b,14);
  a ^= c; a -= ht_rotate_left_u32(c,11);
  b ^= a; b -= ht_rotate_left_u32(a,25);
  c ^= b; c -= ht_rotate_left_u32(b,16);
  a ^= c; a -= ht_rotate_left_u32(c,4);
  b ^= a; b -= ht_rotate_left_u32(a,14);
  c ^= b; c -= ht_rotate_left_u32(b,24);
  return c;
}

static inline uint32_t
ht_finalize2_u32(uint32_t a, uint32_t b)
{
  return ht_finalize3_u32(a, b, 0);
}

static inline uint32_t
gt_uint32_key_mul_hash(uint32_t key)
{
  uint32_t x = ((uint64_t)GOLDEN_RATIO_MULTIPLIER * key) & (~(uint32_t)0);
  return x;
}

static inline uint32_t
gt_uint64_key_mul_hash(uint64_t key)
{
  uint32_t low = key, hi = key >> 32;
  uint32_t x = ht_finalize2_u32(gt_uint32_key_mul_hash(low),
                                gt_uint32_key_mul_hash(hi));
  return x;
}

static inline int
ht_ul_cmp(unsigned long a, unsigned long b)
{
  return (a > b) - (a < b);
}

static inline int
ht_ptr_cmp(void *a, void *b)
{
  return (a > b) - (a < b);
}

#endif
