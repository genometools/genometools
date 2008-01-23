/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/score_matrix.h"
#include "libgtext/luahelper.h"
#include "libgtlua/score_matrix_lua.h"

#define SCOREMATRIX_METATABLE  "GenomeTools.score_matrix"
#define check_score_matrix(L, POS) \
        (ScoreMatrix**) luaL_checkudata(L, POS, SCOREMATRIX_METATABLE)

static int score_matrix_lua_new_read_protein(lua_State *L)
{
  ScoreMatrix **sm;
  const char *path;
  Error *err;
  assert(L);
  path = luaL_checkstring(L, 1);
  sm = lua_newuserdata(L, sizeof (ScoreMatrix*));
  assert(sm);
  err = error_new();
  if (!(*sm = score_matrix_new_read_protein(path, err)))
    return lua_gt_error(L, err); /* handle error */
  error_delete(err);
  assert(*sm);
  luaL_getmetatable(L, SCOREMATRIX_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int score_matrix_lua_get_dimension(lua_State *L)
{
  ScoreMatrix **sm;
  unsigned int dimension;
  sm = check_score_matrix(L, 1);
  dimension = score_matrix_get_dimension(*sm);
  lua_pushinteger(L, dimension);
  return 1;
}

static int score_matrix_lua_get_score(lua_State *L)
{
  ScoreMatrix **sm;
  int idx1, idx2;
  int score;
  sm = check_score_matrix(L, 1);
  idx1 = luaL_checkint(L, 2);
  idx2 = luaL_checkint(L, 3);
  luaL_argcheck(L, idx1 >= 0, 2, "idx1 too small");
  luaL_argcheck(L, idx2 >= 0, 3, "idx2 too small");
  luaL_argcheck(L, idx1 < score_matrix_get_dimension(*sm), 2, "idx1 too large");
  luaL_argcheck(L, idx2 < score_matrix_get_dimension(*sm), 3, "idx2 too large");
  score = score_matrix_get_score(*sm, idx1, idx2);
  lua_pushinteger(L, score);
  return 1;
}

static int score_matrix_lua_delete(lua_State *L)
{
  ScoreMatrix **sm;
  sm = check_score_matrix(L, 1);
  score_matrix_delete(*sm);
  return 0;
}

static const struct luaL_Reg score_matrix_lib_f [] = {
  { "score_matrix_new_read_protein", score_matrix_lua_new_read_protein },
  { NULL, NULL }
};

static const struct luaL_Reg score_matrix_lib_m [] = {
  { "get_dimension", score_matrix_lua_get_dimension },
  { "get_score", score_matrix_lua_get_score },
  { NULL, NULL }
};

int luaopen_score_matrix(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, SCOREMATRIX_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, score_matrix_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, score_matrix_lib_m);
  luaL_register(L, "gt", score_matrix_lib_f);
  return 1;
}
