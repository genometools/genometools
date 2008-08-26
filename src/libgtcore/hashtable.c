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
#include <assert.h>
#include <limits.h>
#include <setjmp.h>
#if TJ_DEBUG > 1
#include <stdio.h>
#endif
#include <string.h>

#include "libgtcore/array.h"
#include "libgtcore/cstr.h"
#include "libgtcore/gt_qsort_r.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"

union link_data
{
  htsize_t *table;
};

#define free_mark (~(htsize_t)0)
#define end_mark (free_mark - 1)
#define mark_bit ((free_mark >> 1) + 1)
enum {
  MIN_SIZE_LOG     =   4,
  FILL_DIVISOR     = 256,
  DEFAULT_LOW_MUL  =  32,       /**< will be used as quotient
                                   DEFAULT_LOW_MUL/FILL_DIVISOR */
  DEFAULT_HIGH_MUL = 192,
};

typedef htsize_t (*GetLinkFunc)(Hashtable *ht, htsize_t idx);
typedef void (*SetLinkFunc)(Hashtable *ht, htsize_t idx, htsize_t link);

static htsize_t
ht_get_table_link(Hashtable *ht, htsize_t idx);
static void
ht_set_table_link(Hashtable *ht, htsize_t idx, htsize_t link);

#if POLYMORPHLINK
#define HT_GET_LINK(ht, idx) ht->get_link(ht, idx)
#define HT_SET_LINK(ht, idx, link) ht->set_link(ht, idx, link)
#else
#define HT_GET_LINK(ht, idx) ht_get_table_link(ht, idx)
#define HT_SET_LINK(ht, idx, link) ht_set_table_link(ht, idx, link)
#endif

struct Hashtable
{
  HashElemInfo table_info;
  void *table;
  GetLinkFunc get_link;
  SetLinkFunc set_link;
  htsize_t table_mask, high_fill, low_fill, current_fill;
  union link_data links;
  unsigned short table_size_log, high_fill_mul, low_fill_mul;
};

static inline void *
ht_elem_ptr(const Hashtable *ht, htsize_t idx)
{
  return (char *)ht->table + ht->table_info.elem_size * idx;
}

static inline void
ht_cp_elem(Hashtable *ht, htsize_t dest_idx, const void *src)
{
  memcpy(ht_elem_ptr(ht, dest_idx), src, ht->table_info.elem_size);
}

static void
ht_resize(Hashtable *ht, unsigned short new_size_log);

extern Hashtable *
hashtable_new(HashElemInfo table_info)
{
  return hashtable_new_with_start_size(table_info, MIN_SIZE_LOG);
}

static void
ht_reinit(Hashtable *ht, HashElemInfo table_info, unsigned short size_log,
        unsigned short high_mul, unsigned short low_mul)
{
  htsize_t table_size;
  assert(high_mul > low_mul);
  assert(low_mul > 0 && high_mul < FILL_DIVISOR);
  ht->table_info = table_info;
  ht->table_size_log = size_log;
  ht->table_mask = (table_size = 1 << size_log) - 1;
  ht->table = ma_realloc(ht->table, table_info.elem_size * table_size);
  ht->high_fill_mul = high_mul;
  ht->high_fill
    = (unsigned long long)ht->high_fill_mul * table_size / FILL_DIVISOR;
  ht->low_fill_mul = low_mul;
  ht->low_fill
    = (unsigned long long)ht->low_fill_mul * table_size / FILL_DIVISOR;
  {
    uint32_t i;
    ht->links.table = ma_realloc(ht->links.table,
                                 sizeof (*(ht->links.table)) * table_size);
    for (i = 0; i < table_size; ++i)
      ht->links.table[i] = free_mark;
    ht->get_link = ht_get_table_link;
    ht->set_link = ht_set_table_link;
  }
}

static void
ht_init(Hashtable *ht, HashElemInfo table_info, unsigned short size_log,
        unsigned short high_mul, unsigned short low_mul)
{
  assert(size_log < sizeof (htsize_t) * CHAR_BIT);
  ht->current_fill = 0;
  ht->table = ht->links.table = NULL;
  ht_reinit(ht, table_info, size_log, high_mul, low_mul);
}

static void
ht_destruct(Hashtable *ht)
{
  ma_free(ht->table);
  ma_free(ht->links.table);
}

extern Hashtable *
hashtable_new_with_start_size(HashElemInfo table_info, unsigned short size_log)
{
  Hashtable *ht;
  ht = ma_malloc(sizeof (*ht));
  ht_init(ht, table_info, size_log, DEFAULT_HIGH_MUL, DEFAULT_LOW_MUL);
  return ht;
}

static int
ht_insert(Hashtable *ht, const void *elem);

static enum iterator_op
ht_insert_wrapper(void *elem, void *data, UNUSED Error *err)
{
#ifndef NDEBUG
  int ins_count =
#endif
    ht_insert(data, elem);
#ifndef NDEBUG
  if (!ins_count)
  {
    fputs("Insertion mysteriously failed on hashtable resize.", stderr);
    abort();
  }
#endif
  return CONTINUE_ITERATION;
}

static void
ht_resize(Hashtable *ht, unsigned short new_size_log)
{
  Hashtable new_ht;
#ifndef NDEBUG
  htsize_t new_size = 1 << new_size_log;
#endif
  assert(ht);
  if (new_size_log != ht->table_size_log)
  {
    ht_init(&new_ht, ht->table_info, new_size_log, ht->high_fill_mul,
            ht->low_fill_mul);
    assert(ht->current_fill < new_size);
    hashtable_foreach(ht, ht_insert_wrapper, &new_ht, NULL);
    ht_destruct(ht);
    memcpy(ht, &new_ht, sizeof (*ht));
  }
}

static inline htsize_t
ht_elem_hash_idx(const Hashtable *ht, const void *elem)
{
  return ht->table_info.keyhash(elem) & ht->table_mask;
}

#define ht_traverse_list_of_key(ht, elem, pre_loop, in_loop, post_loop) \
  do {                                                                  \
    Hashtable *htref = (ht);                                            \
    htsize_t elem_hash = ht_elem_hash_idx(htref, (elem)),               \
      idx, link = elem_hash;                                            \
    pre_loop;                                                           \
    do {                                                                \
      idx = link;                                                       \
      link = HT_GET_LINK(htref, idx);                                   \
      in_loop;                                                          \
    }                                                                   \
    while (!(link & mark_bit));                                         \
    post_loop;                                                          \
  } while (0)

#if TJ_DEBUG > 1
static void
ht_traverse_list_of_key_debug(Hashtable *ht, const void *elem)
{
  Hashtable *htref = (ht);
  htsize_t elem_hash = ht_elem_hash_idx(htref, (elem)),
    idx, link = elem_hash;
  do {
    idx = link;
    link = HT_GET_LINK(htref, idx);
    if (link != free_mark
        && !ht->table_info.cmp(elem, ht_elem_ptr(ht, idx)))
    {
      fputs("found\n", stderr);
    }
  }
  while (!(link & mark_bit));
}
#endif

extern void *
hashtable_get(Hashtable *ht, const void *elem)
{
  assert(ht);
#if TJ_DEBUG > 1
  ht_traverse_list_of_key_debug(ht, elem);
#endif
  ht_traverse_list_of_key(ht, elem, ,
                          if (link != free_mark
                              && !ht->table_info.cmp(elem,
                                                     ht_elem_ptr(ht, idx)))
                            return ht_elem_ptr(ht, idx),);
  return NULL;
}

extern int
hashtable_add(Hashtable *ht, const void *elem)
{
  int insert_count;
  assert(ht && elem);
  if (ht->current_fill + 1 > ht->high_fill)
    ht_resize(ht, ht->table_size_log + 1);
  insert_count = ht_insert(ht, elem);
  return insert_count;
}

static htsize_t
ht_find_free_idx(Hashtable *ht, htsize_t start_idx, int search_dir)
{
  htsize_t new_idx = start_idx;
  assert(ht->current_fill < ht->table_mask + 1);
  do {
    new_idx = (new_idx + search_dir) & ht->table_mask;
  } while (HT_GET_LINK(ht, new_idx) != free_mark);
  return new_idx;
}

static int
ht_insert(Hashtable *ht, const void *elem)
{
  htsize_t insert_pos;
  do {
    htsize_t elem_hash = ht_elem_hash_idx(ht, elem), idx,
      link = elem_hash;
    if (HT_GET_LINK(ht, link) == free_mark)
    {
      /* we can insert at initial link and start a new chain */
      insert_pos = link;
      break;
    }
    else if (ht_elem_hash_idx(ht, ht_elem_ptr(ht, link)) != elem_hash)
    {
      /* relocate chained element */
      htsize_t reloc_idx = link, reloc_referent, new_idx;
      ht_traverse_list_of_key(ht, ht_elem_ptr(ht, reloc_idx),,
                              if (link == reloc_idx)
                                break;,
                              reloc_referent = idx;);
      new_idx = ht_find_free_idx(ht, reloc_referent, -1);
      ht_cp_elem(ht, new_idx, ht_elem_ptr(ht, reloc_idx));
      HT_SET_LINK(ht, new_idx,
                  HT_GET_LINK(ht, reloc_idx));
      HT_SET_LINK(ht, reloc_referent, new_idx);
      insert_pos = link;
      break;
    }
    do {
      idx = link;
      link = HT_GET_LINK(ht, idx);
      if (!ht->table_info.cmp(elem, ht_elem_ptr(ht, idx)))
        /* don't insert elements already present! */
        return 0;
    } while (link != end_mark);
    {
      /* we can search, starting at idx, for a
       * free position to insert elem */
      htsize_t referent = idx,
        new_idx = ht_find_free_idx(ht, idx, +1);
      HT_SET_LINK(ht, referent, new_idx);
      insert_pos = new_idx;
    }
  } while (0);
  ht_cp_elem(ht, insert_pos, elem);
  HT_SET_LINK(ht, insert_pos, end_mark);
  ht->current_fill += 1;
  return 1;
}

static htsize_t
ht_remove(Hashtable *ht, const void *elem);

static inline void
ht_shrink(Hashtable *ht);

extern int
hashtable_remove(Hashtable *ht, const void *elem)
{
  htsize_t remove_pos;
  assert(ht && elem);
  remove_pos = ht_remove(ht, elem);
  if (remove_pos != free_mark)
  {
    ht_shrink(ht);
    return 1;
  }
  else
    return 0;
}

static inline void
ht_shrink(Hashtable *ht)
{
  if (ht->current_fill < ht->low_fill
      && ht->table_size_log > MIN_SIZE_LOG)
  {
    unsigned short new_size_log = ht->table_size_log;
    htsize_t low_fill = ht->low_fill, old_low_fill = low_fill;
    do {
      old_low_fill = low_fill;
      --new_size_log;
      low_fill >>= 1;
    } while (ht->current_fill < old_low_fill && new_size_log > MIN_SIZE_LOG);
    ht_resize(ht, new_size_log);
  }
}

static htsize_t
ht_remove(Hashtable *ht, const void *elem)
{
  htsize_t remove_pos = free_mark, referent = free_mark;
  ht_traverse_list_of_key(ht, elem,,
                          if (link != free_mark
                              && !ht->table_info.cmp(elem,
                                                     ht_elem_ptr(ht, idx)))
                          { remove_pos = idx; break; }
                          referent = idx;,);
  /* was elem found? */
  if (remove_pos != free_mark)
  {
    htsize_t chain_next = HT_GET_LINK(ht, remove_pos);
    if (referent != free_mark)
    {
      HT_SET_LINK(ht, referent, chain_next);
    }
    else if (chain_next != end_mark)
    {
      /* handle removal of chain head, where there is a non-empty tail */
      /* find next free field to move chain head to temporarily */
      htsize_t cp_dest_idx = ht_find_free_idx(ht, remove_pos, -1);
      ht_cp_elem(ht, cp_dest_idx, ht_elem_ptr(ht, remove_pos));
      ht_cp_elem(ht, remove_pos, ht_elem_ptr(ht, chain_next));
      HT_SET_LINK(ht, remove_pos, HT_GET_LINK(ht, chain_next));
      HT_SET_LINK(ht, chain_next, free_mark);
      remove_pos = cp_dest_idx;
    }
    if (ht->table_info.free_op.free_elem_with_data)
      ht->table_info.free_op.free_elem_with_data(ht_elem_ptr(ht, remove_pos),
                                                 ht->table_info.table_data);
    HT_SET_LINK(ht, remove_pos, free_mark);
    --ht->current_fill;
    return remove_pos;
  }
  return free_mark;
}

struct hash_to_array_data
{
  size_t elem_size;
  Array *hash_entries;
};

static enum iterator_op
ht_save_entry_to_array(void *elem, void *data, UNUSED Error *err)
{
  Array *hash_entries;
  assert(elem && data);
  hash_entries = ((struct hash_to_array_data *)data)->hash_entries;
  array_add_elem(hash_entries, elem,
                 ((struct hash_to_array_data *)data)->elem_size);
  return CONTINUE_ITERATION;
}

extern int
hashtable_foreach_ordered(Hashtable *ht, Elemvisitfunc iter, void *data,
                          Compare cmp, Error *err)
{
  Array *hash_entries;
  void *elem;
  unsigned long i;
  int had_err;
  error_check(err);
  assert(ht && iter && cmp);
  hash_entries = array_new(ht->table_info.elem_size);
  {
    struct hash_to_array_data visitor_data = { ht->table_info.elem_size,
                                               hash_entries };
    had_err = hashtable_foreach(ht, ht_save_entry_to_array, &visitor_data,
                                err);
  }
  if (!had_err) {
    size_t hash_size;
    gt_qsort_r(array_get_space(hash_entries), array_size(hash_entries),
               array_elem_size(hash_entries), data, (CompareWithData)cmp);
    hash_size = array_size(hash_entries);
    assert(hash_size == hashtable_fill(ht));
    for (i = 0; !had_err && i < hash_size; i++) {
      elem = array_get(hash_entries, i);
      had_err = iter(elem, data, err);
    }
  }
  array_delete(hash_entries);
  return had_err;
}

extern int
hashtable_foreach_in_default_order(Hashtable *ht, Elemvisitfunc iter,
                                   void *data, Error *err)
{
  return hashtable_foreach_ordered(ht, iter, data, ht->table_info.cmp, err);
}

extern int
hashtable_foreach(Hashtable *ht, Elemvisitfunc visitor, void *data,
                  Error *err)
{
  htsize_t i, table_size = ht->table_mask + 1, deletion_count = 0;
  jmp_buf env;
  while (setjmp(env))
    ;
  for (i = 0; i < table_size; ++i)
  {
    htsize_t idx = i, link = HT_GET_LINK(ht, i);
    /* if i is occupied and forms a chain start */
    if (link != free_mark
        && ht_elem_hash_idx(ht, ht_elem_ptr(ht, i)) == i)
      while (1)
      {
        void *elem;
        elem = ht_elem_ptr(ht, idx);
        switch (visitor(elem, data, err))
        {
        case CONTINUE_ITERATION:
          break;
        case STOP_ITERATION:
          return 0;
        case DELETED_ELEM:
          {
            htsize_t remove_pos = ht_remove(ht, elem);
            ht->table_info.free_op.free_elem_with_data(
              ht_elem_ptr(ht, remove_pos), ht->table_info.table_data);
            ++deletion_count;
          }
          break;
        case MODIFIED_KEY:
          if (ht_elem_hash_idx(ht, elem) != i)
          {
            /* elem now belongs to new chain */
            /* FIXME: handle deferred move */
            fprintf(stderr, "Feature MODIFIED_KEY not implemented yet"
                    " (%s:%d).\n", __FILE__, __LINE__);
            /* push idx */
            abort();
          }
          break;
        case REDO_ITERATION:
          longjmp(env, 1);
          break;
        }
        idx = link;
        if (idx == end_mark)
          break;
        link = HT_GET_LINK(ht, link);
      }
  }
  /* FIXME: re-insert pushed indices and rectify affected chains */
  /* if hashtable shrunk below low_mark, resize */
  if (deletion_count &&
      ht->current_fill < ht->low_fill)
    ht_shrink(ht);
  return 0;
}

extern size_t
hashtable_fill(Hashtable *ht)
{
  assert(ht);
  return ht->current_fill;
}

#define ht_internal_foreach(ht,visitcode)                       \
  do {                                                          \
    htsize_t i, table_size = ht->table_mask + 1;                \
    void *table_data = ht->table_info.table_data;               \
    void *elem = ht->table;                                     \
    size_t elem_size = ht->table_info.elem_size;                \
    if (ht->current_fill)                                       \
      for (i = 0; i < table_size; ++i)                          \
      {                                                         \
        if (HT_GET_LINK(ht, i) != free_mark)                    \
        {                                                       \
          visitcode;                                            \
        }                                                       \
        elem = (char *)elem + elem_size;                        \
      }                                                         \
  } while (0)

extern void
hashtable_reset(Hashtable *ht)
{
  assert(ht);
  FreeFuncWData free_elem_with_data =
    ht->table_info.free_op.free_elem_with_data;
  if (free_elem_with_data)
    ht_internal_foreach(ht, free_elem_with_data(elem, table_data));
  ht->current_fill = 0;
  ht_reinit(ht, ht->table_info, MIN_SIZE_LOG, DEFAULT_HIGH_MUL,
            DEFAULT_LOW_MUL);
}

extern void
hashtable_delete(Hashtable *ht)
{
  if (ht)
  {
    FreeFuncWData free_elem_with_data =
      ht->table_info.free_op.free_elem_with_data;
    if (free_elem_with_data)
      ht_internal_foreach(ht, free_elem_with_data(elem, table_data));
    ht_destruct(ht);
    if (ht->table_info.table_data_free)
      ht->table_info.table_data_free(ht->table_info.table_data);
    ma_free(ht);
  }
}

extern uint32_t
ht_ptr_elem_hash(const void *elem)
{
  /* rotate right by 3 because memory addresses are to often aligned
   * at oct addresses. */
#if CHAR_BIT == 8
  if (sizeof (void *) == 4)
    return uint32_key_mul_hash(ht_rotate_right_u32(*(uint32_t *)elem, 3));
  else if (sizeof (void *) == 8)
    return uint64_key_mul_hash(ht_rotate_right_u64(*(uint64_t *)elem, 3));
#else
#error "pointer size is not a multiple of 8, I'd like to hear of your platform"
#endif
}

extern uint32_t
ht_ul_elem_hash(const void *elem)
{
#if CHAR_BIT == 8
  if (sizeof (void *) == 4)
    return uint32_key_mul_hash(*(uint32_t *)elem);
  else if (sizeof (void *) == 8)
    return uint64_key_mul_hash(*(uint64_t *)elem);
#else
#error "pointer size is not a multiple of 8, I'd like to hear of your platform"
#endif
}

#define ht_u32_mix(a,b,c)                                      \
  {                                                            \
    a -= c;  a ^= ht_rotate_left_u32(c, 4);  c += b;           \
    b -= a;  b ^= ht_rotate_left_u32(a, 6);  a += c;           \
    c -= b;  c ^= ht_rotate_left_u32(b, 8);  b += a;           \
    a -= c;  a ^= ht_rotate_left_u32(c,16);  c += b;           \
    b -= a;  b ^= ht_rotate_left_u32(a,19);  a += c;           \
    c -= b;  c ^= ht_rotate_left_u32(b, 4);  b += a;           \
  }

extern uint32_t
uint32_data_hash(const void *data, size_t length)
{
  uint32_t a,b,c;
  a = b = c = 0xdeadbeef + ((uint32_t)length);
  const uint8_t *k = data;
  /*--------------- all but the last block: affect some 32 bits of (a,b,c) */
  while (length > 12)
  {
    a += k[0];
    a += ((uint32_t)k[1])<<8;
    a += ((uint32_t)k[2])<<16;
    a += ((uint32_t)k[3])<<24;
    b += k[4];
    b += ((uint32_t)k[5])<<8;
    b += ((uint32_t)k[6])<<16;
    b += ((uint32_t)k[7])<<24;
    c += k[8];
    c += ((uint32_t)k[9])<<8;
    c += ((uint32_t)k[10])<<16;
    c += ((uint32_t)k[11])<<24;
    ht_u32_mix(a,b,c);
    length -= 12;
    k += 12;
  }
  /*-------------------------------- last block: affect all 32 bits of (c) */
  switch (length)                   /* all the case statements fall through */
  {
  case 12: c+=((uint32_t)k[11])<<24;
  case 11: c+=((uint32_t)k[10])<<16;
  case 10: c+=((uint32_t)k[9])<<8;
  case 9 : c+=k[8];
  case 8 : b+=((uint32_t)k[7])<<24;
  case 7 : b+=((uint32_t)k[6])<<16;
  case 6 : b+=((uint32_t)k[5])<<8;
  case 5 : b+=k[4];
  case 4 : a+=((uint32_t)k[3])<<24;
  case 3 : a+=((uint32_t)k[2])<<16;
  case 2 : a+=((uint32_t)k[1])<<8;
  case 1 : a+=k[0];
    break;
  case 0 : return c;
  }
  return ht_finalize3_u32(a,b,c);
}

static uint32_t
uint32_str_key_hash(const char *str)
{
  const uint8_t *k = (const uint8_t *)str;
  uint32_t c, h = 0xdeadbeef;
  while ((c = *k++))
    h ^= ((h << 5) + (h >> 2) + c);
  return h;
}

extern uint32_t
ht_cstr_elem_hash(const void *elem)
{
  return uint32_str_key_hash(*(const char **)elem);
}

extern int
ht_ptr_elem_cmp(const void *elemA, const void *elemB)
{
  return ht_ptr_cmp(*(void **)elemA, *(void **)elemB);
}

extern int
ht_ul_elem_cmp(const void *elemA, const void *elemB)
{
  return ht_ul_cmp(*(unsigned long *)elemA, *(unsigned long *)elemB);
}

extern int
ht_cstr_elem_cmp(const void *elemA, const void *elemB)
{
  return strcmp(*(const char **)elemA, *(const char **)elemB);
}

/*
 * link table stuff
 */
static htsize_t
ht_get_table_link(Hashtable *ht, htsize_t idx)
{
  return ht->links.table[idx];
}

static void
ht_set_table_link(Hashtable *ht, htsize_t idx, htsize_t link)
{
  ht->links.table[idx] = link;
}

struct ht_elem_2cstr
{
  char *key, *value;
};

extern void
ht_2ptr_elem_free(void *elem)
{
  struct ht_elem_2cstr *p = elem;
  assert(elem);
  ma_free(p->key);
  ma_free(p->value);
}

#define my_ensure(err_state, predicate)         \
  if (!(predicate)) {                           \
    err_state = -1;                             \
    break;                                      \
  }

static inline void
cstr_cstr_elem_dup(struct ht_elem_2cstr *elem,
                   const char *key_template, const char *value_template)
{
  elem->key = cstr_dup(key_template);
  elem->value = cstr_dup(value_template);
}

static int
hashtable_test(HashElemInfo table_info)
{
  FreeFunc orig_free_elem = table_info.free_op.free_elem;
  char *s1 = "foo", *s2 = "bar";
  Hashtable *ht;
  int had_err = 0;
  struct ht_elem_2cstr elemA = { s1, s2 }, elemB = { s2, s1 };
  table_info.free_op.free_elem = NULL;
  do {
    struct ht_elem_2cstr *elem_p;
    /* empty hash */
    ht = hashtable_new(table_info);
    hashtable_delete(ht);

    /* empty hash with reset */
    ht = hashtable_new(table_info);
    hashtable_reset(ht);
    hashtable_delete(ht);

    /* hashes containing one element */
    ht = hashtable_new(table_info);
    hashtable_add(ht, &elemA);
    my_ensure(had_err, !memcmp(hashtable_get(ht, &elemA), &elemA,
                               table_info.elem_size));
    my_ensure(had_err, !hashtable_get(ht, &elemB));
    hashtable_delete(ht);

    /* hashes containing two elements */
    ht = hashtable_new(table_info);
    hashtable_add(ht, &elemA);
    hashtable_add(ht, &elemB);
    elem_p = hashtable_get(ht, &elemA);
    my_ensure(had_err, elem_p
              && !memcmp(elem_p, &elemA, table_info.elem_size));
    elem_p = hashtable_get(ht, &elemB);
    my_ensure(had_err, elem_p
              && !memcmp(elem_p, &elemB, table_info.elem_size));

    /* remove element A and ensure it's no longer present */
    my_ensure(had_err, hashtable_remove(ht, &elemA));
    my_ensure(had_err, !hashtable_get(ht, &elemA));

    elem_p = hashtable_get(ht, &elemB);
    my_ensure(had_err, elem_p
              && !memcmp(elem_p, &elemB, table_info.elem_size));
    hashtable_delete(ht);

    /* hashes containing two elements (store key and value in
     * hashtable) where simple free is the correct way to get rid of them
     */
    if (orig_free_elem == ht_2ptr_elem_free)
    {
      struct ht_elem_2cstr elem_dup;
      table_info.free_op.free_elem = ht_2ptr_elem_free;
      ht = hashtable_new(table_info);

      cstr_cstr_elem_dup(&elem_dup, s1, s2);
      hashtable_add(ht, &elem_dup);

      cstr_cstr_elem_dup(&elem_dup, s2, s1);
      hashtable_add(ht, &elem_dup);

      elem_p = hashtable_get(ht, &s1);
      my_ensure(had_err, elem_p && !strcmp(elem_p->value, s2));

      elem_p = hashtable_get(ht, &s2);
      my_ensure(had_err, elem_p && !strcmp(elem_p->value, s1));

      my_ensure(had_err, hashtable_remove(ht, &s1)); /* remove first element */
      my_ensure(had_err, !hashtable_get(ht, &s1));

      elem_p = hashtable_get(ht, s2);
      my_ensure(had_err, elem_p && !strcmp(elem_p->value,  s1));
      hashtable_delete(ht);
    }
  } while (0);
  return had_err;
}

int hashtable_unit_test(UNUSED Error *err)
{
  int had_err;
  error_check(err);
  static const HashElemInfo
    hash_ptr = { ht_ptr_elem_hash, { NULL }, sizeof (struct ht_elem_2cstr),
                 ht_ptr_elem_cmp, NULL, NULL },
    hash_str = { ht_cstr_elem_hash, { NULL }, sizeof (struct ht_elem_2cstr),
                 ht_cstr_elem_cmp, NULL, NULL };
  /* hash key as string */
  had_err = hashtable_test(hash_str);

  /* hash key by pointer value */
  if (!had_err)
    had_err = hashtable_test(hash_ptr);

  return had_err;
}
