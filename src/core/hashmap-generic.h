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
#ifndef HASHMAP_GENERIC_H
#define HASHMAP_GENERIC_H

#include <string.h>

#include "core/hashtable.h"
#include "core/unused_api.h"

/**
 * @file hashmap-generic.h
 * @brief Describes a generic set of macros to define
 * - data types
 * - destructor
 * - access functions
 * for a mapping stored in a GtHashtable object
 */

/*
 * In detail, where keytag and valuetag represent the string passed to
 * the macro as the corresponding argument, the following will be
 * declared by DECLARE_HASHMAP:
 *
 * keytag_valuetag_map_entry: struct containing fields of type keytype
 *                            and valuetype respectively
 * keytag_valuetag_gt_hashmap_new: function constructing a hashtable
 *                              setup to store elements of type
 *                              keytag_valuetag_map_entry
 * keytag_valuetag_gt_hashmap_new_with_start_size: function constructing a
 *                            hashtable setup to store elements of type
 *                            keytag_valuetag_map_entry, start_size is
 *                            initially determined
 * keytag_valuetag_gt_hashmap_delete: destructor function corresponding
 *                            to above new
 * keytag_valuetag_gt_hashmap_get: return pointer to value of hashmap
 *                            entry with given key (note: this allows
 *                            to change the hashmap entry in place)
 * keytag_valuetag_gt_hashmap_remove: remove entry with given key from
 *                            hashmap
 * keytag_valuetag_gt_hashmap_add: insert given key and value into the
 *                            hashmap, if no element with key is present
 *                            yet, or update the value of an element
 *                            which compares equal wrt key
 * keytag_valuetag_gt_hashmap_foreach: loops over all entries in hashmap
 *                            and calls visitor function passed by
 *                            function pointer
 * keytag_valuetag_gt_hashmap_foreach_ordered: as _foreach, but in order
 *                            determined by comparison function
 * keytag_valuetag_gt_hashmap_foreach_in_default_order: as _foreach, but
 *                            visit occurs in order defined by
 *                            comparison function used in definition
 *                            of map type
 */

#define NULL_DESTRUCTOR(ptr)

#define DECLARE_HASHMAP(keytype, keytag, valuetype, valuetag, \
                        storagedecl,inlineifstatic)           \
  typedef struct {                                            \
    keytype key;                                              \
    valuetype value;                                          \
  } keytag##_##valuetag##_map_entry;                          \
                                                              \
  storagedecl HashElemInfo keytag##_##valuetag##_hashtype;    \
                                                              \
  typedef int (*keytag##_##valuetag##_gt_hashmap_KeyCmp)(     \
    const keytype a, const keytype b);                        \
                                                              \
  typedef int                                                 \
  (*keytag##_##valuetag##_gt_hashmap_KeyCmpWithData)(         \
    const keytype a, const keytype b, const void *data);      \
                                                              \
  /*@unused@*/ GT_UNUSED static inline                        \
  GtHashtable * keytag##_##valuetag##_gt_hashmap_new(void)    \
  {                                                           \
    return gt_hashtable_new(keytag##_##valuetag##_hashtype);  \
  }                                                           \
                                                              \
  /*@unused@*/ GT_UNUSED static inline                        \
  GtHashtable * keytag##_##valuetag##_gt_hashmap_new_with_start_size(\
                                                unsigned short size_log)  \
  {                                                           \
    return gt_hashtable_new_with_start_size(keytag##_##valuetag##_hashtype,\
                                            size_log);        \
  }                                                           \
                                                              \
  /*@unused@*/ GT_UNUSED static inline void                   \
  keytag##_##valuetag##_gt_hashmap_delete(GtHashtable *ht)    \
  {                                                           \
    gt_hashtable_delete(ht);                                  \
  }                                                           \
                                                              \
  /*@unused@*/ GT_UNUSED static inline valuetype *            \
  keytag##_##valuetag##_gt_hashmap_get(GtHashtable *ht,       \
                                    const keytype key)        \
  {                                                           \
    keytag##_##valuetag##_map_entry *map_entry;               \
    map_entry = gt_hashtable_get(ht, &key);                   \
    if (map_entry != NULL)                                    \
      return &(map_entry->value);                             \
    else                                                      \
      return NULL;                                            \
  }                                                           \
                                                              \
  /*@unused@*/ GT_UNUSED static inline int                    \
  keytag##_##valuetag##_gt_hashmap_remove(GtHashtable *ht,    \
                                       const keytype key)     \
  {                                                           \
    return gt_hashtable_remove(ht, &key);                     \
  }                                                           \
                                                              \
  /*@unused@*/ GT_UNUSED static inline void                   \
  keytag##_##valuetag##_gt_hashmap_add(GtHashtable *ht,       \
                                    const keytype key,        \
                                    valuetype value)          \
  {                                                           \
    void *val;                                                \
    keytag##_##valuetag##_map_entry map_entry                 \
      = { (keytype)key, value };                              \
    if (gt_hashtable_add_with_storage_ptr(ht, &map_entry,     \
                                          &val)) {            \
    } else {                                                  \
      ((keytag##_##valuetag##_map_entry *)                    \
         val)->value = value;                                 \
    }                                                         \
  }                                                           \
                                                              \
  /*@unused@*/ GT_UNUSED static inline void*                  \
  keytag##_##valuetag##_gt_hashmap_add_and_return_storage(GtHashtable *ht,  \
                                    const keytype key,        \
                                    valuetype value)          \
  {                                                           \
    void *val;                                                \
    keytag##_##valuetag##_map_entry map_entry                 \
      = { (keytype)key, value };                              \
    if (gt_hashtable_add_with_storage_ptr(ht, &map_entry,     \
                                          &val)) {            \
    } else {                                                  \
      ((keytag##_##valuetag##_map_entry *)                    \
         val)->value = value;                                 \
    }                                                         \
    return &((keytag##_##valuetag##_map_entry *)              \
              val)->value;                                    \
  }                                                           \
                                                              \
  typedef enum iterator_op                                    \
  (*keytag##_##valuetag##_gt_hashmap_iteratorfunc)(           \
    keytype key, valuetype value, void *data, GtError *err);  \
                                                              \
  /*@unused@*/ GT_UNUSED storagedecl inlineifstatic int       \
  keytag##_##valuetag##_gt_hashmap_foreach(                   \
    GtHashtable *ht,                                          \
    keytag##_##valuetag##_gt_hashmap_iteratorfunc iter,       \
    void *data, GtError *err);                                \
                                                              \
  /*@unused@*/ GT_UNUSED storagedecl inlineifstatic int       \
  keytag##_##valuetag##_gt_hashmap_foreach_ordered(           \
    GtHashtable *ht,                                          \
    keytag##_##valuetag##_gt_hashmap_iteratorfunc iter,       \
    void *data, keytag##_##valuetag##_gt_hashmap_KeyCmp cmp,  \
    GtError *err);                                            \
                                                              \
  /*@unused@*/ GT_UNUSED storagedecl inlineifstatic int       \
  keytag##_##valuetag##_gt_hashmap_foreach_in_default_order(  \
    GtHashtable *ht,                                          \
    keytag##_##valuetag##_gt_hashmap_iteratorfunc iter,       \
    void *data, GtError *err);

#define DECLARE_SAFE_DEREF(valuetype, valuetag)               \
  GT_UNUSED static inline valuetype                           \
  valuetag##_gt_safe_deref(valuetype *ptr)                    \
  {                                                           \
    return ptr ? *ptr : NULL;                                 \
  }

#define DEFINE_HASHMAP(keytype, keytag, valuetype, valuetag,  \
                       keyhash, keycmp, keydestructor,        \
                       valuedestructor, storagedecl,          \
                       inlineifstatic)                        \
                                                              \
  storagedecl inlineifstatic void                             \
  keytag##_##valuetag##_destruct(void *elem)                  \
  {                                                           \
    GT_UNUSED keytag##_##valuetag##_map_entry *map_entry = elem; \
    keydestructor(map_entry->key);                            \
    valuedestructor(map_entry->value);                        \
  }                                                           \
                                                              \
  /*@ignore@*/\
  storagedecl HashElemInfo keytag##_##valuetag##_hashtype = { \
    keyhash,                                                  \
    { keytag##_##valuetag##_destruct },                       \
    sizeof (keytag##_##valuetag##_map_entry),                 \
    keycmp, NULL, NULL                                        \
  };                                                          \
  /*@end@*/\
                                                              \
  typedef struct {                                            \
    void *data;                                               \
    keytag##_##valuetag##_gt_hashmap_iteratorfunc iter;       \
    keytag##_##valuetag##_gt_hashmap_KeyCmp keycmp;           \
  } keytag##_##valuetag##_gt_hashmap_DataIterCmpTripel;       \
                                                              \
  static enum iterator_op                                     \
  keytag##_##valuetag##_gt_hashmap_iter(                      \
    void *elem, void *data, GtError *err)                     \
  {                                                           \
    keytag##_##valuetag##_map_entry *entry = elem;            \
    keytag##_##valuetag##_gt_hashmap_DataIterCmpTripel *dip   \
      = data;                                                 \
    return dip->iter(entry->key, entry->value, dip->data,     \
                     err);                                    \
  }                                                           \
                                                              \
  GT_UNUSED storagedecl inlineifstatic int                    \
  keytag##_##valuetag##_gt_hashmap_foreach(                   \
    GtHashtable *ht,                                          \
    keytag##_##valuetag##_gt_hashmap_iteratorfunc iter,       \
    void *data, GtError *err)                                 \
  {                                                           \
    keytag##_##valuetag##_gt_hashmap_DataIterCmpTripel dip =  \
      { data, iter, NULL };                                   \
    return gt_hashtable_foreach(                              \
      ht, keytag##_##valuetag##_gt_hashmap_iter, &dip, err);  \
  }                                                           \
                                                              \
  GT_UNUSED static int                                        \
  keytag##_##valuetag##_gt_hashmap_cmp(                       \
    const void *elemA, const void *elemB, void *data)         \
  {                                                           \
    const keytag##_##valuetag##_map_entry *entryA = elemA,    \
      *entryB = elemB;                                        \
    keytag##_##valuetag##_gt_hashmap_DataIterCmpTripel *dip   \
      = data;                                                 \
    return ((keytag##_##valuetag##_gt_hashmap_KeyCmpWithData) \
            dip->keycmp)(entryA->key, entryB->key,            \
                         dip->data);                          \
  }                                                           \
                                                              \
  GT_UNUSED storagedecl inlineifstatic int                    \
  keytag##_##valuetag##_gt_hashmap_foreach_ordered(           \
    GtHashtable *ht,                                          \
    keytag##_##valuetag##_gt_hashmap_iteratorfunc iter,       \
    void *data, keytag##_##valuetag##_gt_hashmap_KeyCmp cmp,  \
    GtError *err)                                             \
  {                                                           \
    keytag##_##valuetag##_gt_hashmap_DataIterCmpTripel dip    \
      = { data, iter, cmp };                                  \
    return gt_hashtable_foreach_ordered(                      \
      ht, keytag##_##valuetag##_gt_hashmap_iter, &dip,        \
      (GtCompare)keytag##_##valuetag##_gt_hashmap_cmp, err);  \
  }                                                           \
                                                              \
  GT_UNUSED storagedecl inlineifstatic int                    \
  keytag##_##valuetag##_gt_hashmap_foreach_in_default_order(  \
    GtHashtable *ht,                                          \
    keytag##_##valuetag##_gt_hashmap_iteratorfunc iter,       \
    void *data, GtError *err)                                 \
  {                                                           \
    keytag##_##valuetag##_gt_hashmap_DataIterCmpTripel dip    \
      = { data, iter, NULL };                                 \
    return gt_hashtable_foreach_in_default_order(             \
      ht, keytag##_##valuetag##_gt_hashmap_iter, &dip, err);  \
  }                                                           \

#endif
