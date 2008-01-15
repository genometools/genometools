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

#ifndef LUAHELPER_H
#define LUAHELPER_H

#include "lua.h"
#include "libgtcore/error.h"
#include "libgtcore/strarray.h"

int  lua_set_modules_path(lua_State*, Error*);

void lua_set_arg(lua_State*, const char *argv_0, const char **argv);
void lua_export_metatable(lua_State*, const char *metatable_desc);

void lua_push_strarray_as_table(lua_State*, StrArray*);

/* Propagate the error given in <err> (which must be set) to <L>.
   Takes ownership of the error and deletes it. */
int  lua_gt_error(lua_State *L, Error *err);

#endif
