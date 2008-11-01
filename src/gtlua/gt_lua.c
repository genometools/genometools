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

#include "core/assert_api.h"
#include "gtlua/gtcore_lua.h"
#include "gtlua/gtext_lua.h"
#include "gtlua/gt_lua.h"

#ifndef WITHOUT_CAIRO
#include "gtlua/annotationsketch_lua.h"
#endif

int gt_lua_open_lib(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  gt_lua_open_core(L); /* open core */
  gt_lua_open_extended(L); /* open extended */
#ifndef WITHOUT_CAIRO
  gt_lua_open_annotationsketch(L); /* open annotationsketch */
#endif
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}
