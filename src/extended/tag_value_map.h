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

#ifndef TAG_VALUE_MAP_H
#define TAG_VALUE_MAP_H

#include "core/error.h"

/* A very simple tag/value map absolutely optimized for space (i.e., memory
   consumption) on the cost of time. Basically, each read/write access costs
   O(n) time, whereas n denotes the accumulated length of all tags and values
   contained in the map. Tags and values cannot have length 0.

   The implementation as a char* shines through (also to save one additional
   memory allocation), therefore the usage is a little bit different compared
   to other GenomeTools classes.
   See the implementation of gt_tag_value_map_example() for an example.
*/

typedef char* GtTagValueMap;

typedef void (*GtTagValueMapIteratorFunc)(const char *tag, const char *value,
                                          void *data);

GtTagValueMap gt_tag_value_map_new(const char *tag, const char *value);
void          gt_tag_value_map_delete(GtTagValueMap);
void          gt_tag_value_map_add(GtTagValueMap*, const char *tag,
                                                   const char *value);
const char*   gt_tag_value_map_get(const GtTagValueMap, const char *tag);
void          gt_tag_value_map_foreach(const GtTagValueMap,
                                       GtTagValueMapIteratorFunc,
                                       void *data);
int           gt_tag_value_map_example(GtError*);

#endif
