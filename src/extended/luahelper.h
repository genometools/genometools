/*
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
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
#include "core/error.h"
#include "core/str_array.h"

int   gt_lua_set_modules_path(lua_State*, GtError*);

void  gt_lua_set_arg(lua_State*, const char *argv_0, const char **argv);
void  gt_lua_export_metatable(lua_State*, const char *metatable_desc);

void  gt_lua_push_strarray_as_table(lua_State*, GtStrArray*);
int   gt_lua_get_table_as_strarray(lua_State *L, int index,
                                   GtStrArray *outarray, GtError *err);

/* Propagate the error given in <err> (which must be set) to <L>.
   Takes ownership of the error and deletes it. */
int   gt_lua_error(lua_State *L, GtError *err);
/* Check whether the object at the given stack position <ud> is a userdatum with
   a metatable that matches the given <tname>. It returns NULL if the object
   does not have the correct metatable (or if it is not a userdata); otherwise,
   it returns the userdata address. I.e. like <luaL_checkudata> but without
   throwing an error. */
void* gt_lua_try_checkudata(lua_State *L, int ud, const char *tname);

#endif
