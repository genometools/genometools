/*
  Copyright (c) 2011 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "lauxlib.h"
#include "lua.h"
#include "lualib.h"

#include "match/lua_tools.h"

void gt_lua_stack_dump(lua_State *L) {
  int i;
  int top = lua_gettop(L);
  for (i = 1; i <= top; i++)
  { /* repeat for each level */
    int t = lua_type(L, i);
    switch (t) {

      case LUA_TSTRING:  /* strings */
        fprintf(stderr, "`%s'", lua_tostring(L, i));
        break;

      case LUA_TBOOLEAN:  /* booleans */
        lua_toboolean(L, i) ? fprintf(stderr,"true") : fprintf(stderr,"false");
        break;

      case LUA_TNUMBER:  /* numbers */
        fprintf(stderr, "%g", lua_tonumber(L, i));
        break;

      default:  /* other values */
        fprintf(stderr, "%s", lua_typename(L, t));
        break;
    }
    fprintf(stderr, "  ");  /* put a separator */
  }
  fprintf(stderr, "\n");  /* end the listing */
}
