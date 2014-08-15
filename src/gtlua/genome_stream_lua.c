/*
  Copyright (c) 2007 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  Copyright (c) 2014 Genome Research Ltd

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
#include "core/class_alloc_lock.h"
#include "core/fileutils_api.h"
#include "core/unused_api.h"
#include "extended/node_stream_api.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/luahelper.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/genome_stream_lua.h"

static int gff3_in_stream_lua_new_sorted(lua_State *L)
{
  GtNodeStream **gs;
  const char *filename;
  gt_assert(L);
  /* get/check parameters */
  filename = luaL_checkstring(L, 1);
  luaL_argcheck(L, gt_file_exists(filename), 1, "file does not exist");
  /* construct object */
  gs = lua_newuserdata(L, sizeof (GtNodeStream*));
  *gs = gt_gff3_in_stream_new_sorted(filename);
  gt_assert(*gs);
  luaL_getmetatable(L, GENOME_STREAM_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int gff3_out_stream_lua_new(lua_State *L)
{
  GtNodeStream **out_stream, **in_stream = check_genome_stream(L, 1);
  gt_assert(L);
  /* construct object */
  out_stream = lua_newuserdata(L, sizeof (GtNodeStream*));
  *out_stream = gt_gff3_out_stream_new(*in_stream, NULL);
  gt_assert(*out_stream);
  luaL_getmetatable(L, GENOME_STREAM_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int gff3_out_stream_lua_new_retainids(lua_State *L)
{
  GtNodeStream **out_stream, **in_stream = check_genome_stream(L, 1);
  gt_assert(L);
  /* construct object */
  out_stream = lua_newuserdata(L, sizeof (GtNodeStream*));
  *out_stream = gt_gff3_out_stream_new(*in_stream, NULL);
  gt_assert(*out_stream);
  gt_gff3_out_stream_retain_id_attributes((GtGFF3OutStream*) *out_stream);
  luaL_getmetatable(L, GENOME_STREAM_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int gt_node_stream_lua_next_tree(lua_State *L)
{
  GtNodeStream **gs = check_genome_stream(L, 1);
  GtGenomeNode *gn;
  GtError *err = gt_error_new();
  if (gt_node_stream_next(*gs, &gn, err))
    return gt_lua_error(L, err); /* handle error */
  else if (gn)
    gt_lua_genome_node_push(L, gn);
  else
    lua_pushnil(L);
  gt_error_delete(err);
  return 1;
}

/* This stub only reports an error if trying to run a custom stream with no
   custom behaviour implemented. */
static int gt_node_stream_lua_next_tree_fail(lua_State *L)
{
  GtError *err = gt_error_new();
  gt_error_set(err, "no custom 'next_tree' method defined in custom stream");
  return gt_lua_error(L, err);
}

typedef struct {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  lua_State *L;
  int ref;   /* to allow us to get the stream back when called in the
                middle of a chain */
} GtLuaCustomStream;

const GtNodeStreamClass* gt_lua_custom_stream_class(void);

#define lua_custom_stream_cast(GS)\
        gt_node_stream_cast(gt_lua_custom_stream_class(), GS)

static int lua_custom_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                  GtError *err)
{
  GT_UNUSED GtNodeStream **s;
  GtGenomeNode **retnode;
  GtLuaCustomStream *lcs;
  int had_err = 0;
  gt_assert(ns);
  gt_error_check(err);
  lcs = lua_custom_stream_cast(ns);

  /* push correct Lua object for this stream */
  lua_rawgeti(lcs->L, LUA_REGISTRYINDEX, lcs->ref);
  s = check_genome_stream(lcs->L, -1);

  /* get overridden next_tree method */
  lua_pushstring(lcs->L, "next_tree");
  lua_gettable(lcs->L, -2);

  /* make sure there is one */
  gt_assert(!lua_isnil(lcs->L, -1));
  lua_pushvalue(lcs->L, -2);

  if (lua_pcall(lcs->L, 1, 1, 0)) {
    const char *error = lua_tostring(lcs->L, -1);
    gt_error_set(err, "%s", error);
    had_err = -1;
  }

  if (!had_err) {
    if (lua_isnil(lcs->L, -1)) {
      *gn = NULL;
    } else {
      retnode = gt_lua_try_checkudata(lcs->L, -1, GENOME_NODE_METATABLE);
      if (!retnode) {
        const char *type;
        type = lua_tostring(lcs->L, -1);
        gt_error_set(err, "custom 'next_tree' method must return a genome "
                          "node or nil, was %s", type);
        had_err = -1;
      } else {
        *gn = gt_genome_node_ref(*retnode);
      }
    }
  }
  return had_err;
}

const GtNodeStreamClass* gt_lua_custom_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtLuaCustomStream),
                                   NULL,
                                   lua_custom_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

static int lua_custom_stream_new_generic(lua_State *L, bool sorted)
{
  GtLuaCustomStream *lcs;
  GtNodeStream **ns;
  gt_assert(L);

  ns = lua_newuserdata(L, sizeof (GtNodeVisitor*));
  gt_assert(ns);
  *ns = gt_node_stream_create(gt_lua_custom_stream_class(), sorted);
  gt_assert(*ns);
  lcs = lua_custom_stream_cast(*ns);
  luaL_getmetatable(L, GENOME_STREAM_METATABLE);
  lua_setmetatable(L, -2);

  /* set clean env for this object */
  lua_newtable(L);
  lua_setfenv(L, -2);
  lcs->L = L;

  /* replace the default next_tree method to force an override in the custom
     stream */
  lua_pushstring(L, "next_tree");
  lua_pushcfunction(L, gt_node_stream_lua_next_tree_fail);
  lua_settable(L, -3);

  /* store reference to Lua object */
  lua_pushvalue(L, -1);
  lcs->ref = luaL_ref(L, LUA_REGISTRYINDEX);

  return 1;
}

static int lua_custom_stream_new_sorted(lua_State *L)
{
  return lua_custom_stream_new_generic(L, true);
}

static int lua_custom_stream_new_unsorted(lua_State *L)
{
  return lua_custom_stream_new_generic(L, false);
}

static int gt_node_stream_lua_delete(lua_State *L)
{
  GtNodeStream **gs = check_genome_stream(L, 1);
  gt_node_stream_delete(*gs);
  return 0;
}

static int custom_stream_lua_index(lua_State* L)
{
  lua_getfenv(L, -2);
  lua_pushvalue(L, -2);
  lua_rawget(L, -2);
  if (lua_isnoneornil(L, -1) == 0) {
      return 1;
  }
  lua_pop(L, 2);

  /* check the metatable */
  lua_getmetatable(L, -2);
  lua_pushvalue(L, -2);
  lua_rawget(L, -2);

  return 1;
}

static int custom_stream_lua_newindex(lua_State* L)
{
  lua_getfenv(L, -3);
  lua_pushvalue(L, -3);
  lua_pushvalue(L, -3);
  lua_rawset(L, -3);
  return 0;
}

static const struct luaL_Reg gt_node_stream_lib_f [] = {
  { "gff3_in_stream_new_sorted", gff3_in_stream_lua_new_sorted },
  { "custom_stream_new_sorted", lua_custom_stream_new_sorted },
  { "custom_stream_new_unsorted", lua_custom_stream_new_unsorted },
  { "gff3_out_stream_new", gff3_out_stream_lua_new },
  { "gff3_out_stream_new_retainids", gff3_out_stream_lua_new_retainids },
  { NULL, NULL }
};

static const struct luaL_Reg gt_node_stream_lib_m [] = {
  { "next_tree", gt_node_stream_lua_next_tree },
  { NULL, NULL }
};

int gt_lua_open_genome_stream(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_newmetatable(L, GENOME_STREAM_METATABLE);
  /* metatable.__index = metatable */
  /* lua_pushvalue(L, -1); */ /* duplicate the metatable */
  /* lua_setfield(L, -2, "__index"); */
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, gt_node_stream_lua_delete);
  lua_settable(L, -3);
  lua_pushstring(L, "__index");
  lua_pushcfunction(L, custom_stream_lua_index);
  lua_settable(L, -3);
  lua_pushstring(L, "__newindex");
  lua_pushcfunction(L, custom_stream_lua_newindex);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, gt_node_stream_lib_m);
  lua_pop(L, 1);
  luaL_register(L, "gt", gt_node_stream_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}
