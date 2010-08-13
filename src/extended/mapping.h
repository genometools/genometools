/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef MAPPING_H
#define MAPPING_H

#include "core/str.h"

typedef enum {
  GT_MAPPINGTYPE_STRING,
  GT_MAPPINGTYPE_INTEGER
} GtMappingType;

/* a generic mapping */
typedef struct GtMapping GtMapping ;

/* creates a new mapping from the Lua file <mapping_file> which must contain a
   global table or function with name <global_name>. The global table or
   function must contain entries and return results of the given type,
   respectively. Returns NULL on error. */
GtMapping*  gt_mapping_new(GtStr *mapping_file, const char *global_name,
                           GtMappingType type, GtError*);
/* map <input> to string, returns NULL on error */
GtStr*      gt_mapping_map_string(GtMapping*, const char *input, GtError*);
/* map <input> to integer <output>, returns -1 on error */
int         gt_mapping_map_integer(GtMapping*, long *output, const char *input,
                                   GtError*);
void        gt_mapping_delete(GtMapping*);

#endif
