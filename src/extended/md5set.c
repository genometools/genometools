/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include <ctype.h>
#include "md5.h"
#include "core/log.h"
#include "core/ma.h"
#include "extended/reverse.h"
#include "core/safearith.h"
#include "extended/md5set.h"

typedef struct {
  uint64_t l, h;
} gt_md5_t;

#define GT_MD5_T_EQUAL(A,B) \
  ((A).l == (B).l && (A).h == (B).h)

#define GT_MD5_T_IS_EMPTY(A) \
  ((A).l == 0 && (A).h == 0)

static const unsigned long long gt_md5set_primes[] =
{
#include "extended/md5set.prtab"
};

#define GT_MD5SET_TOO_LARGE \
  "fatal: no prime number larger than %llu in lookup table\n" \
  "developers: modify scripts/makeprimestable.sh to create a larger table\n"

/* if n is in lookup table return it;
 * otherwise return first element in lookup table larger than n */
static unsigned long gt_md5set_get_size(unsigned long n)
{
  unsigned long long u, l, i, k, k_i;

  k = (unsigned long long)n;
  if (k > GT_MD5SET_LARGEST_PRIME)
  {
    fprintf(stderr, GT_MD5SET_TOO_LARGE, k);
    exit(1);
  }
  if (k < gt_md5set_primes[0])
    return (unsigned long)gt_md5set_primes[0];
  l = 0;
  u = GT_MD5SET_NOFPRIMES - 1;
  do {
    i = (l + u) >> 1;
    gt_assert(u >= l);
    if (u - l == 1)
    {
      return (gt_md5set_primes[l] == k)
          ? (unsigned long)k
          : (unsigned long)gt_md5set_primes[u];
    }
    k_i = gt_md5set_primes[i];
    if (k < k_i)
    {
      u = i;
    }
    else if (k > k_i)
    {
      l = i;
    }
    else
    {
      return (unsigned long)k_i;
    }
  } while (1);
}

struct GtMd5set
{
  /* hash table */
  gt_md5_t       *table;
  unsigned long  alloc;
  unsigned long  fill;
  unsigned long  maxfill;
  /* string temp buffer */
  char           *buffer;
  unsigned long  bufsize;
};

#define GT_MD5SET_PREPARE_INSERTION(SET) \
  ((SET)->fill)++;\
  if ((SET)->fill > (SET)->maxfill) \
    gt_md5set_alloc_table((SET), \
        /* using alloc + 1, the next value in the lookup table is returned */ \
        gt_md5set_get_size((SET)->alloc + 1));\
  gt_assert((SET)->fill <= (SET)->maxfill)

static void gt_md5set_alloc_table(GtMd5set *set, unsigned long newsize);

enum GtMd5setSearchResult
{
  GT_MD5SET_EMPTY,
  GT_MD5SET_KEY_FOUND,
  GT_MD5SET_COLLISION
};

static inline enum GtMd5setSearchResult gt_md5set_search_pos(GtMd5set *set,
    gt_md5_t k, bool insert_if_not_found, unsigned long i)
{
  if (GT_MD5_T_IS_EMPTY(set->table[i]))
  {
    if (insert_if_not_found)
    {
      GT_MD5SET_PREPARE_INSERTION(set);
      set->table[i] = k;
    }
    return GT_MD5SET_EMPTY;
  }
  return (GT_MD5_T_EQUAL(k, set->table[i]))
         ? GT_MD5SET_KEY_FOUND
         : GT_MD5SET_COLLISION;
}

#define GT_MD5SET_H1(MD5, TABLE_SIZE) \
  ((MD5).l % (TABLE_SIZE))

/* h2 result must be relatively prime to table_size;
   any value 1 < x <= table_size - 1 is it, because table_size is prime */
#define GT_MD5SET_H2(MD5, TABLE_SIZE) \
  (((MD5).h % ((TABLE_SIZE) - 1)) + 1)

static bool gt_md5set_search(GtMd5set *set, gt_md5_t k,
    bool insert_if_not_found)
{
  unsigned long i, c;
#ifndef NDEBUG
  unsigned long first_i;
#endif
  enum GtMd5setSearchResult retval;

  i = GT_MD5SET_H1(k, set->alloc);
  retval = gt_md5set_search_pos(set, k, insert_if_not_found, i);
  if (retval != GT_MD5SET_COLLISION)
    return retval;

  /* open addressing by double hashing */
  c = GT_MD5SET_H2(k, set->alloc);
  gt_assert(c > 0);
#ifndef NDEBUG
  first_i = i;
#endif
  while (1)
  {
    i = (i + c) % set->alloc;
    gt_assert(i != first_i);
    retval = gt_md5set_search_pos(set, k, insert_if_not_found, i);
    if (retval != GT_MD5SET_COLLISION)
      return retval;
  }
}

static void gt_md5set_rehash(GtMd5set *set, gt_md5_t *oldtable,
    unsigned long oldsize)
{
  unsigned long i;
  set->fill = 0;
  for (i = 0; i < oldsize; i++)
    if (!GT_MD5_T_IS_EMPTY(oldtable[i]))
      gt_md5set_search(set, oldtable[i], true);
}

#define GT_MD5SET_MAX_LOAD_FACTOR 0.8

static void gt_md5set_alloc_table(GtMd5set *set, unsigned long newsize)
{
  gt_md5_t *oldtable;
  unsigned long oldsize;

  oldsize = set->alloc;
  oldtable = set->table;
  set->alloc = newsize;
  set->maxfill =
    (unsigned long)((double)newsize * GT_MD5SET_MAX_LOAD_FACTOR);
  set->table = gt_calloc(newsize, sizeof (gt_md5_t));
  if (oldtable != NULL)
  {
    gt_log_log("rehashing %lu elements; old size: %lu, new size: %lu\n",
        set->fill, oldsize, newsize);
    gt_md5set_rehash(set, oldtable, oldsize);
    gt_free(oldtable);
  }
}

GtMd5set *gt_md5set_new(unsigned long number_of_elements)
{
  GtMd5set *md5set;
  md5set = gt_malloc(sizeof (GtMd5set));
  md5set->fill = 0;
  md5set->alloc = 0;
  md5set->table = NULL;
  gt_md5set_alloc_table(md5set, gt_md5set_get_size(number_of_elements +
        (number_of_elements >> 2)));
  gt_assert(number_of_elements < md5set->maxfill);
  md5set->buffer = NULL;
  md5set->bufsize = 0;
  return md5set;
}

void gt_md5set_delete(GtMd5set *md5set)
{
  if (md5set != NULL)
  {
    gt_free(md5set->table);
    gt_free(md5set->buffer);
    gt_free(md5set);
  }
}

static void gt_md5set_prepare_buffer(GtMd5set *md5set, unsigned long bufsize)
{
  gt_assert(md5set != NULL);

  if (md5set->buffer == NULL)
  {
    md5set->buffer = gt_malloc(sizeof (char) * bufsize);
    md5set->bufsize = bufsize;
  }
  else if (md5set->bufsize < bufsize)
  {
    md5set->buffer = gt_realloc(md5set->buffer, sizeof (char) * bufsize);
    md5set->bufsize = bufsize;
  }
}

#define GT_MD5SET_NOT_FOUND  0
#define GT_MD5SET_FOUND      1
#define GT_MD5SET_RC_FOUND   2

#define GT_MD5SET_HASH_STRING(BUF, LEN, MD5) \
  md5((BUF), gt_safe_cast2long(LEN), (char*)&(MD5))

int gt_md5set_add_sequence(GtMd5set *md5set, const char* seq,
    unsigned long seqlen, bool double_strand, GtError *err)
{
  gt_md5_t md5sum, md5sum_rc;
  unsigned long i;
  int retval = 0;
  bool found;

  gt_assert(md5set != NULL);
  gt_assert(md5set->table != NULL);

  gt_md5set_prepare_buffer(md5set, seqlen);
  for (i = 0; i < seqlen; i++)
    md5set->buffer[i] = toupper(seq[i]);

  GT_MD5SET_HASH_STRING(md5set->buffer, seqlen, md5sum);
  found = gt_md5set_search(md5set, md5sum, true);
  if (found)
    return GT_MD5SET_FOUND;

  if (double_strand)
  {
    retval = gt_reverse_complement(md5set->buffer, seqlen, err);
    if (retval != 0)
    {
      gt_assert(retval < 0);
      return retval;
    }

    GT_MD5SET_HASH_STRING(md5set->buffer, seqlen, md5sum_rc);
    found = gt_md5set_search(md5set, md5sum_rc, false);
    if (found)
      return GT_MD5SET_RC_FOUND;
  }

  return GT_MD5SET_NOT_FOUND;
}
