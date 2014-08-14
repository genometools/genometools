/*
  Copyright (c) 2008-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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

#include <stdlib.h>
#include <string.h>
#include "core/ma.h"
#include "core/ensure.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/tag_value_map.h"

/* The GtTagValueMap is implemented as a simple char* which points to a memory
   region organized as follows:

   tag\0value\0tag\0value\0\0
*/

GtTagValueMap gt_tag_value_map_new(const char *tag, const char *value)
{
  GtTagValueMap map;
  size_t tag_len, value_len;
  gt_assert(tag && value);
  tag_len = strlen(tag);
  value_len = strlen(value);
  gt_assert(tag_len && value_len);
  map = gt_malloc((tag_len + 1 + value_len + 1 + 1) * sizeof *map);
  memcpy(map, tag, tag_len + 1);
  memcpy(map + tag_len + 1, value, value_len + 1);
  map[tag_len + 1 + value_len + 1] = '\0';
  return map;
}

/* Stores map length in <map_len> if the return value equals NULL (i.e., if not
   value has been found) and <map_len> does not equal NULL. */
static char* get_value(const GtTagValueMap map, const char *tag,
                       size_t *map_len)
{
  char *map_ptr;
  const char *tag_ptr;
  /* search for equal tag */
  map_ptr = map;
  tag_ptr = tag;
  for (;;) {
    while (*map_ptr == *tag_ptr++) {
      if (*map_ptr++ == '\0')
        return map_ptr; /* match found -> return value */
    }
    /* no match found */
    while (*map_ptr++ != '\0'); /* go to next value */
    if (*map_ptr == '\0')
      break; /* no next value found */
    /* reset tag pointer */
    tag_ptr = tag;
  }
  if (map_len)
    *map_len = map_ptr - map; /* save map length */
  return NULL;
}

static size_t get_map_len(const GtTagValueMap map)
{
  const char *map_ptr = map;
  for (;;) {
    if ((*map_ptr++ == '\0') && (*map_ptr++ == '\0'))
      break;
  }
  return map_ptr - map - 1;
}

static size_t get_map_nof_items(const GtTagValueMap map)
{
  const char *map_ptr = map;
  GtUword nof_items = 0;
  bool seen = false;
  for (;;) {
    if (*map_ptr++ == '\0') {
      if (*(map_ptr) == '\0') {
        if (seen) {
          nof_items++;
        }
        break;
      }
      if (seen) {
        nof_items++;
        seen = false;
      } else {
        seen = true;
      }
    }
  }
  return nof_items;
}

void gt_tag_value_map_add(GtTagValueMap *map, const char *tag,
                          const char *value)
{
  size_t tag_len, value_len, map_len = 0;
  GT_UNUSED const char *tag_already_used;
  gt_assert(map && *map && tag && value);
  tag_len = strlen(tag);
  value_len = strlen(value);
  gt_assert(tag_len && value_len);
  /* determine current map length */
  tag_already_used = get_value(*map, tag, &map_len);
  gt_assert(!tag_already_used); /* map does not contain given <tag> already */
  /* allocate additional space */
  *map = gt_realloc(*map, map_len + tag_len + 1 + value_len + 1 + 1);
  /* store new tag/value pair */
  memcpy(*map + map_len, tag, tag_len + 1);
  memcpy(*map + map_len + tag_len + 1, value, value_len + 1);
  (*map)[map_len + tag_len + 1 + value_len + 1] = '\0';
}

void gt_tag_value_map_remove(GtTagValueMap *map, const char *tag)
{
  size_t tag_len, value_len, map_len;
  char *value;
  gt_assert(map && tag && get_map_nof_items(*map) > 1);
  tag_len = strlen(tag);
  gt_assert(tag_len);
  /* determine value for given tag */
  value = get_value(*map, tag, NULL);
  gt_assert(value);
  map_len = get_map_len(*map);
  value_len = strlen(value);
  /* move memory from end position of value to start position of tag */
  memmove(value - tag_len - 1, value + value_len + 1,
          map_len - ((size_t) value - (size_t) *map + value_len));
  *map = gt_realloc(*map, map_len - (tag_len + 1 + value_len + 1) + 1);
  gt_assert((*map)[map_len - (tag_len + 1 + value_len + 1)] == '\0');
}

void gt_tag_value_map_set(GtTagValueMap *map, const char *tag,
                          const char *new_value)
{
  size_t old_value_len, new_value_len, map_len = 0;
  char *old_value;
  gt_assert(map && *map && tag && new_value);
  gt_assert(strlen(tag));
  new_value_len = strlen(new_value);
  gt_assert(new_value_len);
  /* determine current map length */
  old_value = get_value(*map, tag, &map_len);
  if (!old_value)
    return gt_tag_value_map_add(map, tag, new_value);
  /* tag already used -> replace it */
  old_value_len = strlen(old_value);
  map_len = get_map_len(*map);
  if (new_value_len < old_value_len) {
    memcpy(old_value, new_value, new_value_len);
    memmove(old_value + new_value_len, old_value + old_value_len,
            map_len - ((size_t) old_value - (size_t) *map + old_value_len) + 1);
    *map = gt_realloc(*map, map_len - (old_value_len - new_value_len) + 1);
  }
  else if (new_value_len == old_value_len) {
    memcpy(old_value, new_value, new_value_len);
  }
  else { /* (new_value_len > old_value_len)  */
    *map = gt_realloc(*map, map_len + (new_value_len - old_value_len) + 1);
    /* determine old_value again, realloc() might have moved it */
    old_value = get_value(*map, tag, &map_len);
    gt_assert(old_value);
    memmove(old_value + new_value_len, old_value + old_value_len,
            map_len - ((size_t) old_value - (size_t) *map + old_value_len) + 1);
    memcpy(old_value, new_value, new_value_len);
  }
  gt_assert((*map)[map_len - old_value_len + new_value_len] == '\0');
}

const char* gt_tag_value_map_get(const GtTagValueMap map, const char *tag)
{
  gt_assert(map && tag && strlen(tag));
  return get_value(map, tag, NULL);
}

GtUword gt_tag_value_map_size(const GtTagValueMap map)
{
  gt_assert(map);
  return get_map_nof_items(map);
}

void gt_tag_value_map_foreach(const GtTagValueMap map,
                              GtTagValueMapIteratorFunc func,
                              void *data)
{
  const char *map_ptr, *tag;
  gt_assert(map && func);
  map_ptr = map;
  do { /* the map has at least one tag/value pair */
    tag = map_ptr;
    while (*map_ptr++ != '\0'); /* skip tag and \0 */
    func(tag, map_ptr /* value */, data);
    while (*map_ptr++ != '\0'); /* skip value and \0 */
  } while (*map_ptr != '\0');
}

void gt_tag_value_map_show(const GtTagValueMap map)
{
  bool null_terminator_read = false;
  const char *map_ptr;
  gt_assert(map);
  map_ptr = map;
  for (;;) {
    if (*map_ptr == '\0') {
      printf("\\0");
      if (null_terminator_read)
        break;
      else
        null_terminator_read = true;
    }
    else {
      gt_xputchar(*map_ptr);
      null_terminator_read = false;
    }
    map_ptr++;
  }
  gt_xputchar('\n');
}

int gt_tag_value_map_example(GT_UNUSED GtError *err)
{
  GtTagValueMap map;

  gt_error_check(err);

  map = gt_tag_value_map_new("tag 1", "value 1");
  gt_tag_value_map_add(&map, "tag 2", "value 2");
  gt_tag_value_map_add(&map, "tag 3", "value 3");

  gt_assert(!gt_tag_value_map_get(map, "unused tag"));
  gt_assert(!strcmp(gt_tag_value_map_get(map, "tag 1"), "value 1"));
  gt_assert(!strcmp(gt_tag_value_map_get(map, "tag 2"), "value 2"));
  gt_assert(!strcmp(gt_tag_value_map_get(map, "tag 3"), "value 3"));

  gt_tag_value_map_delete(map);

  return 0;
}

static GtTagValueMap create_filled_tag_value_list(void)
{
  GtTagValueMap map = gt_tag_value_map_new("tag 1", "value 1");
  gt_tag_value_map_add(&map, "tag 2", "value 2");
  gt_tag_value_map_add(&map, "tag 3", "value 3");

  gt_assert(!gt_tag_value_map_get(map, "unused tag"));
  gt_assert(!strcmp(gt_tag_value_map_get(map, "tag 1"), "value 1"));
  gt_assert(!strcmp(gt_tag_value_map_get(map, "tag 2"), "value 2"));
  gt_assert(!strcmp(gt_tag_value_map_get(map, "tag 3"), "value 3"));

  return map;
}

int gt_tag_value_map_unit_test(GtError *err)
{
  GtTagValueMap map;
  int had_err = 0;

  gt_error_check(err);

  /* test gt_tag_value_map_set() (new tags are shorter than old tags) */
  map = create_filled_tag_value_list();
  gt_ensure(gt_tag_value_map_size(map) == 3);
  gt_tag_value_map_set(&map, "tag 1", "val X");
  gt_tag_value_map_set(&map, "tag 2", "val Y");
  gt_tag_value_map_set(&map, "tag 3", "val Z");
  gt_ensure(!gt_tag_value_map_get(map, "unused tag"));
  gt_ensure(!strcmp(gt_tag_value_map_get(map, "tag 1"), "val X"));
  gt_ensure(!strcmp(gt_tag_value_map_get(map, "tag 2"), "val Y"));
  gt_ensure(!strcmp(gt_tag_value_map_get(map, "tag 3"), "val Z"));
  gt_tag_value_map_delete(map);

  /* test gt_tag_value_map_set() (new tags have same length) */
   if (!had_err) {
    map = create_filled_tag_value_list();
    gt_tag_value_map_set(&map, "tag 1", "value X");
    gt_tag_value_map_set(&map, "tag 2", "value Y");
    gt_tag_value_map_set(&map, "tag 3", "value Z");
    gt_ensure(!gt_tag_value_map_get(map, "unused tag"));
    gt_ensure(!strcmp(gt_tag_value_map_get(map, "tag 1"), "value X"));
    gt_ensure(!strcmp(gt_tag_value_map_get(map, "tag 2"), "value Y"));
    gt_ensure(!strcmp(gt_tag_value_map_get(map, "tag 3"), "value Z"));
    gt_tag_value_map_delete(map);
  }

  /* test gt_tag_value_map_set() (new tags are longer than old tags) */
  if (!had_err) {
    map = create_filled_tag_value_list();
    gt_tag_value_map_set(&map, "tag 1", "value XXX");
    gt_tag_value_map_set(&map, "tag 2", "value YYY");
    gt_tag_value_map_set(&map, "tag 3", "value ZZZ");
    gt_ensure(!gt_tag_value_map_get(map, "unused tag"));
    gt_ensure(
              !strcmp(gt_tag_value_map_get(map, "tag 1"), "value XXX"));
    gt_ensure(
              !strcmp(gt_tag_value_map_get(map, "tag 2"), "value YYY"));
    gt_ensure(
              !strcmp(gt_tag_value_map_get(map, "tag 3"), "value ZZZ"));
    gt_tag_value_map_delete(map);
  }

  /* test gt_tag_value_map_remove() (remove first tag)*/
  if (!had_err) {
    size_t old_map_len, new_map_len;
    map = create_filled_tag_value_list();
    gt_tag_value_map_set(&map, "tag 1", "value XXX");
    gt_tag_value_map_set(&map, "tag 2", "value YYY");
    gt_tag_value_map_set(&map, "tag 3", "value ZZZ");
    old_map_len = get_map_len(map);
    gt_tag_value_map_remove(&map, "tag 1");
    gt_ensure(gt_tag_value_map_size(map) == 2);
    new_map_len = get_map_len(map);
    gt_ensure(!gt_tag_value_map_get(map, "unused tag"));
    gt_ensure(!(old_map_len - new_map_len - strlen("tag 1") - 1
                     - strlen("value XXX") -1));
    gt_tag_value_map_delete(map);
  }

  /* test gt_tag_value_map_remove() (remove middle tag)*/
  if (!had_err) {
    size_t old_map_len, new_map_len;
    map = create_filled_tag_value_list();
    gt_tag_value_map_set(&map, "tag 1", "value XXX");
    gt_tag_value_map_set(&map, "tag 2", "value YYY");
    gt_tag_value_map_set(&map, "tag 3", "value ZZZ");
    old_map_len = get_map_len(map);
    gt_tag_value_map_remove(&map, "tag 2");
    gt_ensure(gt_tag_value_map_size(map) == 2);
    new_map_len = get_map_len(map);
    gt_ensure(!gt_tag_value_map_get(map, "unused tag"));
    gt_ensure(!(old_map_len - new_map_len - strlen("tag 2") - 1
                     - strlen("value YYY") -1));
    gt_tag_value_map_delete(map);
  }

  /* test gt_tag_value_map_remove() (remove last tag)*/
  if (!had_err) {
    size_t old_map_len, new_map_len;
    map = create_filled_tag_value_list();
    gt_tag_value_map_set(&map, "tag 1", "value XXX");
    gt_tag_value_map_set(&map, "tag 2", "value YYY");
    gt_tag_value_map_set(&map, "tag 3", "value ZZZ");
    old_map_len = get_map_len(map);
    gt_tag_value_map_remove(&map, "tag 3");
    gt_ensure(gt_tag_value_map_size(map) == 2);
    new_map_len = get_map_len(map);
    gt_ensure(!gt_tag_value_map_get(map, "unused tag"));
    gt_ensure(!(old_map_len - new_map_len - strlen("tag 3") - 1
                     - strlen("value ZZZ") -1));
    gt_tag_value_map_delete(map);
  }

  return had_err;
}

void gt_tag_value_map_delete(GtTagValueMap map)
{
  if (!map) return;
  gt_free(map);
}
