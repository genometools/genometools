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

#include "lauxlib.h"
#include "core/error.h"
#include "extended/luahelper.h"
#include "extended/stream_evaluator.h"
#include "gtlua/genome_stream_lua.h"
#include "gtlua/genome_visitor_lua.h"
#include "gtlua/gtext_lua.h"
#include "gtlua/stream_evaluator_lua.h"

#define STREAM_EVALUATOR_METATABLE  "GenomeTools.stream_evaluator"
#define check_stream_evaluator(L) \
        (GtStreamEvaluator**) luaL_checkudata(L, 1, STREAM_EVALUATOR_METATABLE)

static int stream_evaluator_lua_new(lua_State *L)
{
  GtStreamEvaluator **stream_evaluator;
  GtNodeStream **reality_stream, **prediction_stream;
  reality_stream = check_genome_stream(L, 1);
  prediction_stream = check_genome_stream(L, 2);
  stream_evaluator = lua_newuserdata(L, sizeof (GtStreamEvaluator*));
  *stream_evaluator = gt_stream_evaluator_new(*reality_stream,
                                              *prediction_stream,
                                              true, false, 0);
  luaL_getmetatable(L, STREAM_EVALUATOR_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int stream_evaluator_lua_evaluate(lua_State *L)
{
  GtStreamEvaluator **stream_evaluator;
  GtNodeVisitor **genome_visitor;
  GtError *err;
  stream_evaluator = check_stream_evaluator(L);
  if (lua_gettop(L) >= 2) {
    genome_visitor = check_genome_visitor(L, 2);
  }
  else
    genome_visitor = NULL;
  err = gt_error_new();
  if (gt_stream_evaluator_evaluate(*stream_evaluator, false, false, false,
                                   genome_visitor ? *genome_visitor : NULL,
                                   err)) {
    return gt_lua_error(L, err);
  }
  gt_error_delete(err);
  return 0;
}

static int stream_evaluator_lua_show(lua_State *L)
{
  GtStreamEvaluator **stream_evaluator = check_stream_evaluator(L);
  gt_stream_evaluator_show(*stream_evaluator, NULL);
  return 0;
}

static int stream_evaluator_lua_delete(lua_State *L)
{
  GtStreamEvaluator **stream_evaluator;
  stream_evaluator = check_stream_evaluator(L);
  gt_stream_evaluator_delete(*stream_evaluator);
  return 0;
}

static const struct luaL_Reg stream_evaluator_lib_f [] = {
  { "stream_evaluator_new", stream_evaluator_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg stream_evaluator_lib_m [] = {
  { "evaluate", stream_evaluator_lua_evaluate },
  { "show", stream_evaluator_lua_show },
  { NULL, NULL }
};

int gt_lua_open_stream_evaluator(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_newmetatable(L, STREAM_EVALUATOR_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, stream_evaluator_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, stream_evaluator_lib_m);
  lua_pop(L, 1);
  luaL_register(L, "gt", stream_evaluator_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}
