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

/*@unused@*/ static inline uint32_t
gt_ht_rotate_left_u32(uint32_t m, unsigned short k)
{
  return (((m)<<(k)) | ((m)>>((sizeof (m) * CHAR_BIT)-(k))));
}

/*@unused@*/ static inline uint64_t
gt_ht_rotate_left_u64(uint64_t m, unsigned short k)
{
  return (((m)<<(k)) | ((m)>>((sizeof (m) * CHAR_BIT)-(k))));
}

/*@unused@*/ static inline uint32_t
gt_ht_rotate_riggt_ht_u32(uint32_t m, unsigned short k)
{
  return (((m)>>(k)) | ((m)<<((sizeof (m) * CHAR_BIT)-(k))));
}

/*@unused@*/ static inline uint64_t
gt_ht_rotate_riggt_ht_u64(uint64_t m, unsigned short k)
{
  return (((m)>>(k)) | ((m)<<((sizeof (m) * CHAR_BIT)-(k))));
}

/*
 * The finalize function and the multiplicative hash functions is
 * derived from Bob Jenkins' public domain implementation
 * (http://www.burtleburtle.net/bob/)
 */

static inline uint32_t
gt_ht_finalize3_u32(uint32_t a, uint32_t b, uint32_t c)
{
  c ^= b; c -= gt_ht_rotate_left_u32(b,(unsigned short) 14);
  a ^= c; a -= gt_ht_rotate_left_u32(c,(unsigned short) 11);
  b ^= a; b -= gt_ht_rotate_left_u32(a,(unsigned short) 25);
  c ^= b; c -= gt_ht_rotate_left_u32(b,(unsigned short) 16);
  a ^= c; a -= gt_ht_rotate_left_u32(c,(unsigned short) 4);
  b ^= a; b -= gt_ht_rotate_left_u32(a,(unsigned short) 14);
  c ^= b; c -= gt_ht_rotate_left_u32(b,(unsigned short) 24);
  return c;
}

/*@unused@*/ static inline uint32_t
gt_ht_finalize2_u32(uint32_t a, uint32_t b)
{
  return gt_ht_finalize3_u32(a, b, 0);
}

/*@unused@*/ static inline uint32_t
gt_uint32_key_mul_hash(uint32_t key)
{
  uint32_t x = ((uint64_t)GOLDEN_RATIO_MULTIPLIER * key) & (~(uint32_t)0);
  return x;
}

/*@unused@*/ static inline uint32_t
gt_uint64_key_mul_hash(uint64_t key)
{
  uint32_t low = key, hi = key >> 32;
  uint32_t x = gt_ht_finalize2_u32(gt_uint32_key_mul_hash(low),
                                gt_uint32_key_mul_hash(hi));
  return x;
}

static inline int
gt_ht_ul_cmp(unsigned long a, unsigned long b)
{
  return (int) (a > b) - (int) (a < b);
}

/*@unused@*/ static inline int
gt_ht_ptr_cmp(void *a, void *b)
{
  return (int) (a > b) - (int) (a < b);
}

#endif
