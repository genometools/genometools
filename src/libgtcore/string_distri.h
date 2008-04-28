/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef STRING_DISTRI_H
#define STRING_DISTRI_H

#include "libgtcore/error.h"
#include "libgtcore/genfile.h"

/* A discrete distribution */
typedef struct StringDistri StringDistri;

typedef void (*StringDistriIterFunc)(const char *string,
                                     unsigned long occurrences,
                                     double probability, void *data);

StringDistri*      string_distri_new(void);
void               string_distri_add(StringDistri*, const char*);
/* <sd> must contain at least on element with given <key>. */
void               string_distri_sub(StringDistri *sd, const char *key);
unsigned long      string_distri_get(const StringDistri*, const char*);
void               string_distri_foreach(const StringDistri*,
                                        StringDistriIterFunc, void *data);
void               string_distri_delete(StringDistri*);

#endif
