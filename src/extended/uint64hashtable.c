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

#include "core/assert_api.h"
#include "core/ensure.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/hashtable.h" /* gt_uint64_key_mul_hash */
#include "core/spacecalc.h"
#include "extended/uint64hashtable.h"
#include "extended/uint64hashtable_primes.h"

#define GT_UINT64TABLE_TOO_LARGE \
  "fatal: no prime number larger than %lu in lookup table\n" \
  "developers: modify scripts/makeprimestable.sh to create a larger table\n"

#define GT_UINT64TABLE_NOFPRIMES \
  (sizeof (gt_uint64hashtable_primes) / sizeof (gt_uint64hashtable_primes[0]))

#define GT_UINT64TABLE_LARGEST_PRIME \
  gt_uint64hashtable_primes[GT_UINT64TABLE_NOFPRIMES - 1UL]

/* if n is in lookup table return it;
 * otherwise return first element in lookup table larger than n */
static size_t gt_uint64hashtable_get_size(size_t n)
{
  size_t u, l, i, k, k_i;

  k = n;
  if (k > (size_t)GT_UINT64TABLE_LARGEST_PRIME)
  {
    fprintf(stderr, GT_UINT64TABLE_TOO_LARGE, (unsigned long)k);
    exit(1);
  }
  if (k < (size_t)gt_uint64hashtable_primes[0])
  {
    return (size_t)gt_uint64hashtable_primes[0];
  }
  l = 0;
  u = GT_UINT64TABLE_NOFPRIMES - (size_t)1;
  do {
    i = (l + u) >> 1;
    gt_assert(u >= l);
    if (u - l == (size_t)1)
    {
      return ((size_t)gt_uint64hashtable_primes[l] == k)
          ? k
          : (size_t)gt_uint64hashtable_primes[u];
    }
    k_i = (size_t)gt_uint64hashtable_primes[i];
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
      return k_i;
    }
  } while (true);
  gt_assert(false); /* this point should never be reached */
  return 0;
}

typedef struct
{
  uint64_t key;
  unsigned long count;
} GtUint64hashstoredvalue;

struct GtUint64hashtable
{
  GtUint64hashstoredvalue *hspace;
  size_t    alloc, fill, maxfill;
  unsigned long countcollision;
  bool      zero_occurs;
};

static void gt_uint64hashtable_alloc_table(GtUint64hashtable *table,
                                           size_t tsize);

#define GT_UINT64TABLE_MAX_LOAD_FACTOR 0.8

GtUint64hashtable *gt_uint64hashtable_new(size_t nof_elements)
{
  GtUint64hashtable *table;
  table = gt_malloc(sizeof (*table));
  table->fill = 0;
  table->alloc = 0;
  table->countcollision = 0;
  table->hspace = NULL;
  table->zero_occurs = false;
  gt_uint64hashtable_alloc_table(table, gt_uint64hashtable_get_size(
        (size_t)(1 + (double)nof_elements / GT_UINT64TABLE_MAX_LOAD_FACTOR)));
  gt_assert(nof_elements < table->maxfill);
  return table;
}

void gt_uint64hashtable_delete(GtUint64hashtable *table)
{
  if (table != NULL)
  {
    printf("# number of collisions %lu\n",table->countcollision);
    gt_free(table->hspace);
    gt_free(table);
  }
}

#define GT_UINT64TABLE_EMPTYMARK 0UL

enum GtUint64hashtableSearchResult
{
  GT_UINT64TABLE_EMPTY,
  GT_UINT64TABLE_KEY_FOUND,
  GT_UINT64TABLE_COLLISION
};

static inline enum GtUint64hashtableSearchResult
           gt_uint64hashtable_search_pos(GtUint64hashtable *table,
                                         uint64_t key,
                                         bool insert_if_not_found,
                                         size_t pos)
{
  gt_assert(pos < table->alloc);
  if (table->hspace[pos].key == GT_UINT64TABLE_EMPTYMARK)
  {
    if (insert_if_not_found)
    {
      table->fill++;
      if (table->fill > table->maxfill)
      {
        /* using alloc + 1, the next value in the lookup table is returned */
        gt_uint64hashtable_alloc_table(
                    table,
                    gt_uint64hashtable_get_size(table->alloc + 1UL));
      }
      gt_assert(table->fill <= table->maxfill);
      table->hspace[pos].key = key;
      table->hspace[pos].count = 1UL;
    }
    return GT_UINT64TABLE_EMPTY;
  } else
  {
    if (table->hspace[pos].key == key)
    {
      table->hspace[pos].count++;
      return GT_UINT64TABLE_KEY_FOUND;
    }
    return GT_UINT64TABLE_COLLISION;
  }
}

static inline size_t gt_uint64hashtable_h1(uint64_t key, size_t table_size)
{
  return (size_t) gt_uint64_key_mul_hash(key) % table_size;
}

static inline size_t gt_uint64hashtable_h2(uint64_t key, size_t table_size)
{
  return (size_t) 1 +
         (size_t) gt_uint64_key_mul_hash(key) % (table_size - 1);
}

bool gt_uint64hashtable_search(GtUint64hashtable *table, uint64_t key,
                               bool insert_if_not_found)
{
  gt_assert(table != NULL);
  if (key == 0)
  {
    if (table->zero_occurs)
    {
      return true;
    } else
    {
      if (insert_if_not_found)
      {
        table->zero_occurs = true;
      }
      return false;
    }
  } else
  {
    size_t i, c;
#ifndef NDEBUG
    size_t first_i;
#endif
    enum GtUint64hashtableSearchResult retval;

    i = gt_uint64hashtable_h1(key, table->alloc);
    retval = gt_uint64hashtable_search_pos(table, key, insert_if_not_found, i);
    if (retval != GT_UINT64TABLE_COLLISION)
    {
      return retval == GT_UINT64TABLE_EMPTY ? false : true;
    }
    table->countcollision++;
    c = gt_uint64hashtable_h2(key, table->alloc);
    gt_assert(c > 0);
#ifndef NDEBUG
    first_i = i;
#endif
    while (1)
    {
      i += c;
      if (i > table->alloc)
      {
        i -= table->alloc;
      }
      gt_assert(i != first_i);
      retval = gt_uint64hashtable_search_pos(table, key,insert_if_not_found,i);
      if (retval != GT_UINT64TABLE_COLLISION)
      {
        return retval == GT_UINT64TABLE_EMPTY ? false : true;
      }
    }
  }
}

static void gt_uint64hashtable_rehash(GtUint64hashtable *table,
                                      GtUint64hashstoredvalue *oldtable,
                                      size_t oldsize)
{
  size_t i;

  table->fill = 0;
  for (i = 0; i < oldsize; i++)
  {
    if (oldtable[i].key != GT_UINT64TABLE_EMPTYMARK)
    {
      (void) gt_uint64hashtable_search(table, oldtable[i].key, true);
    }
  }
}

static void gt_uint64hashtable_alloc_table(GtUint64hashtable *table,
                                           size_t tsize)
{
  GtUint64hashstoredvalue *oldspace;
  size_t oldsize;

  oldsize = table->alloc;
  oldspace = table->hspace;
  table->alloc = tsize;
  table->maxfill = (size_t)((double)tsize * GT_UINT64TABLE_MAX_LOAD_FACTOR);
  table->hspace = gt_calloc(tsize, sizeof (*table->hspace));
  printf("calloc (%.2f MB)\n",GT_MEGABYTES(tsize * sizeof (*table->hspace)));
  if (oldspace != NULL)
  {
    gt_log_log("rehashing %lu elements; old size: %lu, new size: %lu\n",
        (unsigned long) table->fill, (unsigned long) oldsize,
        (unsigned long) tsize);
    gt_uint64hashtable_rehash(table, oldspace, tsize);
    gt_free(oldspace);
  }
}

int gt_uint64hashtable_unit_test(GtError *err)
{
  int had_err = 0;
  GtUint64hashtable *table = NULL;
  bool found;
  size_t i, nof_elements;

  gt_error_check(err);

  table = gt_uint64hashtable_new(0);
  ensure(had_err, table != NULL);
  found = gt_uint64hashtable_search(table, (uint64_t)7, false);
  ensure(had_err, !found);
  found = gt_uint64hashtable_search(table, (uint64_t)7, true);
  ensure(had_err, !found);
  found = gt_uint64hashtable_search(table, (uint64_t)7, true);
  ensure(had_err, found);
  gt_uint64hashtable_delete(table);

  nof_elements = (size_t)10000;
  table = gt_uint64hashtable_new(nof_elements);
  ensure(had_err, table != NULL);
  for (i = 0; i < nof_elements; i++)
  {
    found = gt_uint64hashtable_search(table, (uint64_t)i, true);
    ensure(had_err, !found);
  }
  for (i = 0; i < nof_elements; i++)
  {
    found = gt_uint64hashtable_search(table, (uint64_t)i, true);
    ensure(had_err, found);
  }
  gt_uint64hashtable_delete(table);
  return had_err;
}
