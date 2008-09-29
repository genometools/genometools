/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <stdlib.h>
#include <string.h>
#include "core/ma.h"
#include "core/unused_api.h"
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

void gt_tag_value_map_delete(GtTagValueMap map)
{
  if (!map) return;
  gt_free(map);
}

/* Stores map length in <map_len> if the return value equals NULL (i.e., if not
   value has been found) and <map_len> does not equal NULL. */
static const char* get_value(const GtTagValueMap map, const char *tag,
                             size_t *map_len)
{
  const char *map_ptr, *tag_ptr;
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

void gt_tag_value_map_add(GtTagValueMap *map, const char *tag,
                          const char *value)
{
  size_t tag_len, value_len, map_len = 0;
  const char *tag_already_used;
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

const char* gt_tag_value_map_get(const GtTagValueMap map, const char *tag)
{
  gt_assert(map && tag && strlen(tag));
  return get_value(map, tag, NULL);
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
