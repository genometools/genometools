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

#include "libgtcore/str.h"

typedef enum {
  MAPPINGTYPE_STRING,
  MAPPINGTYPE_INTEGER
} MappingType;

/* a generic mapping */
typedef struct Mapping Mapping;

/* creates a new mapping from the Lua file <mapping_file> which must contain a
   global table or function with name <global_name>. The global table or
   function must contain entries and return results of the given type,
   respectively. Returns NULL on error. */
Mapping* mapping_new(Str *mapping_file, const char *global_name,
                     MappingType type, Error*);
/* map <input> to string, returns NULL on error */
Str*     mapping_map_string(Mapping*, const char *input, Error*);
/* map <input> to integer <output>, returns -1 on error */
int      mapping_map_integer(Mapping*, long *output, const char *input, Error*);
void     mapping_delete(Mapping*);

#endif
