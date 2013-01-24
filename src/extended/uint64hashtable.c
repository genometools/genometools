/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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
#include "core/format64.h"
#include "core/qsort_r_api.h"
#include "core/timer_api.h"
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

static int compareGtUint64hashstoredvalue(const void *a,const void *b,
                                          void *data)
{
  GtUint64hashstoredvalue *hspace = (GtUint64hashstoredvalue *) data;
  uint32_t va = *(const uint32_t *) a;
  uint32_t vb = *(const uint32_t *) b;

  if (hspace[va].key < hspace[vb].key)
  {
    return -1;
  }
  if (hspace[va].key > hspace[vb].key)
  {
    return 1;
  }
  gt_assert(false);
  return 0;
}

struct GtUint64hashtable
{
  GtUint64hashstoredvalue *hspace;
  uint32_t *sortedhspace;
  size_t alloc;
  unsigned long countcollisions, zero_count,
                allentries;
  bool zero_occurs;
};

GtUint64hashtable *gt_uint64hashtable_new(size_t nof_elements)
{
  GtUint64hashtable *table;
  const double loadfactor = 1.30;

  table = gt_malloc(sizeof (*table));
  table->countcollisions = 0;
  table->allentries = 0;
  table->zero_occurs = false;
  table->zero_count = 0;
  table->alloc
    = gt_uint64hashtable_get_size((size_t)(1 + loadfactor *
                                               (double) nof_elements));
  table->hspace = gt_calloc(table->alloc, sizeof (*table->hspace));
  return table;
}

void gt_uint64hashtable_delete(GtUint64hashtable *table)
{
  if (table != NULL)
  {
    gt_free(table->hspace);
    gt_free(table);
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
  if (key > 0)
  {
    size_t pos, hashadd = 0, iteration;
    const uint64_t emptymark = 0;

#ifndef NDEBUG
    size_t first_pos;
#endif
    pos = gt_uint64hashtable_h1(key, table->alloc);
#ifndef NDEBUG
    first_pos = pos;
#endif
    for (iteration = 0; iteration < table->alloc; iteration++)
    {
      gt_assert(pos < table->alloc);
      if (table->hspace[pos].key == emptymark)
      {
        if (insert_if_not_found)
        {
          table->allentries++;
          table->hspace[pos].key = key;
          table->hspace[pos].count++;
        }
        return false;
      }
      if (table->hspace[pos].key == key)
      {
        gt_assert(table->hspace[pos].count > 0);
        table->hspace[pos].count++;
        return true;
      }
      table->countcollisions++;
      if (hashadd == 0)
      {
        hashadd = gt_uint64hashtable_h2(key, table->alloc);
      }
      gt_assert(hashadd > 0);
      pos += hashadd;
      if (pos >= table->alloc)
      {
        pos -= table->alloc;
      }
      gt_assert(pos != first_pos);
    }
    fprintf(stderr, "function %s, file %s, line %d.\n"
                    "Cannot find empty slot in hashtable: "
                    "This is probably a bug, please report it.\n",
                    __func__, __FILE__, __LINE__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  } else
  {
    if (!table->zero_occurs)
    {
      if (insert_if_not_found)
      {
        table->zero_occurs = true;
        table->zero_count++;
      }
      return false;
    } else
    {
      gt_assert(table->zero_count > 0);
      table->zero_count++;
      return true;
    }
  }
}

unsigned long gt_uint64hashtable_insertionindex(GtUint64hashtable *table,
                                                uint64_t key)
{
  gt_assert(table != NULL);
  if (key > 0)
  {
    size_t pos, hashadd = 0, iteration;
    const uint64_t emptymark = 0;

#ifndef NDEBUG
    size_t first_pos;
#endif
    pos = gt_uint64hashtable_h1(key, table->alloc);
#ifndef NDEBUG
    first_pos = pos;
#endif
    for (iteration = 0; iteration < table->alloc; iteration++)
    {
      gt_assert(pos < table->alloc);
      if (table->hspace[pos].key == emptymark)
      {
        return ULONG_MAX;
      }
      if (table->hspace[pos].key == key)
      {
        gt_assert(table->hspace[pos].count > 0);
        return --table->hspace[pos].count;
      }
      table->countcollisions++;
      if (hashadd == 0)
      {
        hashadd = gt_uint64hashtable_h2(key, table->alloc);
      }
      gt_assert(hashadd > 0);
      pos += hashadd;
      if (pos >= table->alloc)
      {
        pos -= table->alloc;
      }
      gt_assert(pos != first_pos);
    }
    fprintf(stderr, "function %s, file %s, line %d.\n"
                    "Cannot find empty slot in hashtable: "
                    "This is probably a bug, please report it.\n",
                    __func__, __FILE__, __LINE__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  } else
  {
    gt_assert(table->zero_occurs);
    gt_assert(table->zero_count > 0);
    return --table->zero_count;
  }
}

unsigned long gt_uint64hashtable_countsum_get(const GtUint64hashtable *table)
{
  size_t idx;
  unsigned long sumcount = 0;

  for (idx=0; idx < table->alloc; idx++)
  {
    if (table->hspace[idx].count > 0)
    {
      sumcount += table->hspace[idx].count;
    }
  }
  return sumcount + table->zero_count;
}

unsigned long gt_uint64hashtable_partialsums(GtUint64hashtable *table,
                                             GtTimer *timer)
{
  size_t idx, next = 0;
  unsigned long psum, maxsize = 0;

  table->sortedhspace = gt_malloc((size_t) table->allentries *
                                  sizeof (*table->sortedhspace));
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "sorting the hashkeys",stdout);
  }
  for (idx = 0; idx < table->alloc; idx++)
  {
    if (table->hspace[idx].count > 0)
    {
      gt_assert(next < (size_t) table->allentries);
      table->sortedhspace[next++] = idx;
      if (maxsize < table->hspace[idx].count)
      {
        maxsize = table->hspace[idx].count;
      }
    }
  }
  gt_qsort_r(table->sortedhspace,next,sizeof (*table->sortedhspace),
             table->hspace,compareGtUint64hashstoredvalue);
  gt_assert(next > 0);
  if (table->zero_occurs)
  {
    table->hspace[table->sortedhspace[0]].count += table->zero_count;
  }
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "computing partial sums",stdout);
  }
  for (idx = (size_t) 1; idx < next; idx++)
  {
    table->hspace[table->sortedhspace[idx]].count +=
      table->hspace[table->sortedhspace[idx-1]].count;
  }
  psum = table->hspace[table->sortedhspace[next-1]].count;
  gt_free(table->sortedhspace);
  return psum;
}

int gt_uint64hashtable_unit_test(GtError *err)
{
  int had_err = 0;
  GtUint64hashtable *table = NULL;
  bool found;
  size_t i, nof_elements;

  gt_error_check(err);

  table = gt_uint64hashtable_new(0);
  gt_ensure(had_err, table != NULL);
  found = gt_uint64hashtable_search(table, (uint64_t)7, false);
  gt_ensure(had_err, !found);
  found = gt_uint64hashtable_search(table, (uint64_t)7, true);
  gt_ensure(had_err, !found);
  found = gt_uint64hashtable_search(table, (uint64_t)7, true);
  gt_ensure(had_err, found);
  gt_uint64hashtable_delete(table);

  nof_elements = (size_t)10000;
  table = gt_uint64hashtable_new(nof_elements);
  gt_ensure(had_err, table != NULL);
  for (i = 0; i < nof_elements; i++)
  {
    found = gt_uint64hashtable_search(table, (uint64_t)i, true);
    gt_ensure(had_err, !found);
  }
  for (i = 0; i < nof_elements; i++)
  {
    found = gt_uint64hashtable_search(table, (uint64_t)i, true);
    gt_ensure(had_err, found);
  }
  gt_uint64hashtable_delete(table);
  return had_err;
}
