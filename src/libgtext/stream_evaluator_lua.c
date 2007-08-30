/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "lauxlib.h"
#include "gtlua.h"
#include "libgtcore/env.h"
#include "libgtext/genome_stream_lua.h"
#include "libgtext/genome_visitor_lua.h"
#include "libgtext/gtext_lua.h"
#include "libgtext/stream_evaluator.h"
#include "libgtext/stream_evaluator_lua.h"

#define STREAM_EVALUATOR_METATABLE  "GenomeTools.stream_evaluator"
#define check_stream_evaluator(L) \
        (StreamEvaluator**) luaL_checkudata(L, 1, STREAM_EVALUATOR_METATABLE);

static int stream_evaluator_lua_new(lua_State *L)
{
  StreamEvaluator **stream_evaluator;
  GenomeStream **reality_stream, **prediction_stream;
  Env *env = get_env_from_registry(L);
  reality_stream = check_genome_stream(L, 1);
  prediction_stream = check_genome_stream(L, 2);
  stream_evaluator = lua_newuserdata(L, sizeof (StreamEvaluator**));
  *stream_evaluator = stream_evaluator_new(*reality_stream, *prediction_stream,
                                           false, 0, env);
  luaL_getmetatable(L, STREAM_EVALUATOR_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int stream_evaluator_lua_evaluate(lua_State *L)
{
  StreamEvaluator **stream_evaluator;
  GenomeVisitor **genome_visitor;
  Env *env = get_env_from_registry(L);
  stream_evaluator = check_stream_evaluator(L);
  if (lua_gettop(L) >= 2) {
    genome_visitor = check_genome_visitor(L, 2);
  }
  else
    genome_visitor = NULL;
  stream_evaluator_evaluate(*stream_evaluator, false, false,
                            genome_visitor ?  *genome_visitor : NULL, env);
  return 0;
}

static int stream_evaluator_lua_show(lua_State *L)
{
  StreamEvaluator **stream_evaluator = check_stream_evaluator(L);
  stream_evaluator_show(*stream_evaluator, stdout);
  return 0;
}

static int stream_evaluator_lua_delete(lua_State *L)
{
  StreamEvaluator **stream_evaluator;
  Env *env = get_env_from_registry(L);
  stream_evaluator = check_stream_evaluator(L);
  stream_evaluator_delete(*stream_evaluator, env);
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

int luaopen_stream_evaluator(lua_State *L)
{
  assert(L);
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
  luaL_register(L, "gt", stream_evaluator_lib_f);
  return 1;
}
