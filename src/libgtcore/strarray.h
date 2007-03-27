/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef STRARRAY_H
#define STRARRAY_H

#include <libgtcore/env.h>

/* the string array class  */
typedef struct StrArray StrArray;

StrArray*     strarray_new(Env*);
void          strarray_add_cstr(StrArray*, const char*, Env*);
const char*   strarray_get(const StrArray*, unsigned long strnum);
unsigned long strarray_size(const StrArray*); /* returns number of strings */
void          strarray_delete(StrArray*, Env*);

#endif
