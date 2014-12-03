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
#include "core/compat.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/safearith.h"
#include "core/types_api.h"
#include "extended/md5set.h"
#include "extended/md5set_primes_table.h"
#include "extended/reverse_api.h"

typedef struct {
  uint64_t l, h;
} md5_t;

#define MD5_T_EQUAL(A, B) \
        ((A).l == (B).l && (A).h == (B).h)

#define MD5_T_IS_EMPTY(A) \
        ((A).l == 0 && (A).h == 0)

#define MD5SET_TOO_LARGE \
  "fatal: no prime number larger than "GT_LLU" in lookup table\n" \
  "developers: modify scripts/makeprimestable.sh to create a larger table\n"

#ifndef S_SPLINT_S
/* if n is in lookup table return it;
 * otherwise return first element in lookup table larger than n */
static GtUword md5set_get_size(GtUword n)
{
  GtUint64 u, l, i, k, k_i;

  k = (GtUint64)n;
  if (k > GT_MD5SET_LARGEST_PRIME) {
    fprintf(stderr, MD5SET_TOO_LARGE, k);
    exit(1);
  }
  if (k < gt_md5set_primes[0])
    return (GtUword)gt_md5set_primes[0];
  l = 0;
  u = GT_MD5SET_NOFPRIMES - 1;
  do {
    i = (l + u) >> 1;
    gt_assert(u >= l);
    if (u - l == 1) {
      return (gt_md5set_primes[l] == k)
          ? (GtUword)k
          : (GtUword)gt_md5set_primes[u];
    }
    k_i = gt_md5set_primes[i];
    if (k < k_i) {
      u = i;
    }
    else if (k > k_i) {
      l = i;
    }
    else {
      return (GtUword)k_i;
    }
  } while (1);
  return 0;
}
#endif /* S_SPLINT_S */

struct GtMD5Set
{
  /* hash table */
  md5_t    *table;
  GtUword  alloc;
  GtUword  fill;
  GtUword  maxfill;
  /* string temp buffer */
  char           *buffer;
  GtUword  bufsize;
};

#define GT_MD5SET_PREPARE_INSERTION(SET) \
  ((SET)->fill)++;\
  if ((SET)->fill > (SET)->maxfill) \
    md5set_alloc_table((SET), \
        /* using alloc + 1, the next value in the lookup table is returned */ \
        md5set_get_size((SET)->alloc + 1));\
  gt_assert((SET)->fill <= (SET)->maxfill)

static void md5set_alloc_table(GtMD5Set *set, GtUword newsize);

enum MD5SetSearchResult
{
  GT_MD5SET_EMPTY,
  GT_MD5SET_KEY_FOUND,
  GT_MD5SET_COLLISION
};

static inline enum MD5SetSearchResult md5set_search_pos(GtMD5Set *set,
                                                        md5_t k,
                                                        bool
                                                        insert_if_not_found,
                                                        GtUword i)
{
  if (MD5_T_IS_EMPTY(set->table[i])) {
    if (insert_if_not_found) {
      GT_MD5SET_PREPARE_INSERTION(set);
      set->table[i] = k;
    }
    return GT_MD5SET_EMPTY;
  }
  return (MD5_T_EQUAL(k, set->table[i]))
         ? GT_MD5SET_KEY_FOUND
         : GT_MD5SET_COLLISION;
}

#define MD5SET_H1(MD5, TABLE_SIZE) \
        ((MD5).l % (TABLE_SIZE))

/* h2 result must be relatively prime to table_size;
   any value 1 < x <= table_size - 1 is it, because table_size is prime */
#define MD5SET_H2(MD5, TABLE_SIZE) \
        (((MD5).h % ((TABLE_SIZE) - 1)) + 1)

static bool md5set_search(GtMD5Set *set, md5_t k, bool insert_if_not_found)
{
  GtUword i, c;
#ifndef NDEBUG
  GtUword first_i;
#endif
  enum MD5SetSearchResult retval;

  i = (GtUword) MD5SET_H1(k, set->alloc);
  retval = md5set_search_pos(set, k, insert_if_not_found, i);
  if (retval != GT_MD5SET_COLLISION)
    return retval == GT_MD5SET_EMPTY ? false : true;

  /* open addressing by double hashing */
  c = (GtUword) MD5SET_H2(k, set->alloc);
  gt_assert(c > 0);
#ifndef NDEBUG
  first_i = i;
#endif
  while (1) {
    i = (i + c) % set->alloc;
    gt_assert(i != first_i);
    retval = md5set_search_pos(set, k, insert_if_not_found, i);
    if (retval != GT_MD5SET_COLLISION)
      return retval == GT_MD5SET_EMPTY ? false : true;
  }
}

static void md5set_rehash(GtMD5Set *set, md5_t *oldtable, GtUword oldsize)
{
  GtUword i;
  set->fill = 0;
  for (i = 0; i < oldsize; i++)
    if (!MD5_T_IS_EMPTY(oldtable[i]))
      (void) md5set_search(set, oldtable[i], true);
}

#define MD5SET_MAX_LOAD_FACTOR 0.8

static void md5set_alloc_table(GtMD5Set *set, GtUword newsize)
{
  md5_t *oldtable;
  GtUword oldsize;

  oldsize = set->alloc;
  oldtable = set->table;
  set->alloc = newsize;
  set->maxfill =
    (GtUword)((double)newsize * MD5SET_MAX_LOAD_FACTOR);
  set->table = gt_calloc((size_t)newsize, sizeof (md5_t));
  if (oldtable != NULL) {
    gt_log_log("rehashing " GT_WU " elements; old size: " GT_WU ", new size: "
               GT_WU "\n", set->fill, oldsize, newsize);
    md5set_rehash(set, oldtable, oldsize);
    gt_free(oldtable);
  }
}

GtMD5Set *gt_md5set_new(GtUword nof_elements)
{
  GtMD5Set *md5set;
  md5set = gt_malloc(sizeof (GtMD5Set));
  md5set->fill = 0;
  md5set->alloc = 0;
  md5set->table = NULL;
  md5set_alloc_table(md5set,
                     md5set_get_size(nof_elements + (nof_elements >> 2)));
  gt_assert(nof_elements < md5set->maxfill);
  md5set->buffer = NULL;
  md5set->bufsize = 0;
  return md5set;
}

void gt_md5set_delete(GtMD5Set *set)
{
  if (set != NULL) {
    gt_free(set->table);
    gt_free(set->buffer);
    gt_free(set);
  }
}

static void md5set_prepare_buffer(GtMD5Set *md5set, GtUword bufsize)
{
  gt_assert(md5set != NULL);

  if (md5set->buffer == NULL) {
    md5set->buffer = gt_malloc(sizeof (char) * bufsize);
    md5set->bufsize = bufsize;
  }
  else if (md5set->bufsize < bufsize) {
    md5set->buffer = gt_realloc(md5set->buffer, sizeof (char) * bufsize);
    md5set->bufsize = bufsize;
  }
}

#define MD5SET_HASH_STRING(BUF, LEN, MD5) \
        md5((BUF), gt_safe_cast2long(LEN), (char*)&(MD5))

GtMD5SetStatus gt_md5set_add_sequence(GtMD5Set *set, const char* seq,
                                      GtUword seqlen, bool both_strands,
                                      GtError *err)
{
  md5_t md5sum, md5sum_rc;
  GtUword i;
  int retval = 0;
  bool found;

  gt_assert(set != NULL);
  gt_assert(set->table != NULL);

  md5set_prepare_buffer(set, seqlen);
  for (i = 0; i < seqlen; i++)
    set->buffer[i] = toupper(seq[i]);

  MD5SET_HASH_STRING(set->buffer, seqlen, md5sum);
  found = md5set_search(set, md5sum, true);
  if (found)
    return GT_MD5SET_FOUND;

  if (both_strands) {
    retval = gt_reverse_complement(set->buffer, seqlen, err);
    if (retval != 0) {
      gt_assert(retval < 0);
      return GT_MD5SET_ERROR;
    }

    MD5SET_HASH_STRING(set->buffer, seqlen, md5sum_rc);
    /* if the MD5 sum of the reverse complement equals the MD5 sum of the
       sequence itself we don't check if the reverse complement is in the set.
       Otherwise such sequences would never be added to the set at all. */
    if (md5sum_rc.l == md5sum.l && md5sum_rc.h == md5sum.h) {
      return GT_MD5SET_NOT_FOUND;
    }
    found = md5set_search(set, md5sum_rc, false);
    if (found)
      return GT_MD5SET_RC_FOUND;
  }

  return GT_MD5SET_NOT_FOUND;
}
