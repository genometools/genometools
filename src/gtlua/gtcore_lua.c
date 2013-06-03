/*
  Copyright (c) 2007-2009 Gordon Gremme <gordon@gremme.org>
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
#include "gtlua/alphabet_lua.h"
#include "gtlua/bittab_lua.h"
#include "gtlua/encseq_lua.h"
#include "gtlua/gtcore_lua.h"
#include "gtlua/mathsupport_lua.h"
#include "gtlua/range_lua.h"
#include "gtlua/score_matrix_lua.h"
#include "gtlua/translate_lua.h"

int gt_lua_open_core(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  gt_lua_open_alphabet(L);
  gt_lua_open_bittab(L);
  gt_lua_open_encseq(L);
  gt_lua_open_mathsupport(L);
  gt_lua_open_range(L);
  gt_lua_open_score_matrix(L);
  gt_lua_open_translate(L);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}
