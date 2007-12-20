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

#ifndef RANGE_LUA_H
#define RANGE_LUA_H

#include "lua.h"
#include "libgtcore/range.h"

/* exports the Range class to Lua:

   range       = gt.range_new(startpos, endpos)
   startpos    = range:get_start()
   endpos      = range:get_end()
   range_table = gt.ranges_sort(range_table)
   boolean     = gt.ranges_are_sorted(range_table)
*/
int luaopen_range(lua_State*);

/* push a Range to Lua, returns 1 (number of elements pushed) */
int range_lua_push(lua_State*, Range);

#define RANGE_METATABLE  "GenomeTools.range"
#define check_range(L, POS) \
          (Range*) luaL_checkudata(L, POS, RANGE_METATABLE)

#endif
