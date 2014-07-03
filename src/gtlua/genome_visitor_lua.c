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
#include "core/error_api.h"
#include "core/unused_api.h"
#include "extended/feature_node_api.h"
#include "extended/gff3_visitor_api.h"
#include "extended/node_visitor_api.h"
#include "extended/luahelper.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/genome_visitor_lua.h"

typedef struct {
  const GtNodeVisitor parent_instance;
  lua_State *L;
} GtLuaCustomVisitor;

const GtNodeVisitorClass* lua_custom_visitor_class(void);

#define lua_custom_visitor_cast(GV)\
        gt_node_visitor_cast(lua_custom_visitor_class(), GV)

static int gff3_visitor_lua_new(lua_State *L)
{
  GtNodeVisitor **gv;
  gt_assert(L);
  /* construct object */
  gv = lua_newuserdata(L, sizeof (GtNodeVisitor*));
  *gv = gt_gff3_visitor_new(NULL);
  gt_assert(*gv);
  luaL_getmetatable(L, GENOME_VISITOR_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int lua_custom_visitor_visit_node_generic(GT_UNUSED GtNodeVisitor *nv,
                                                 GtGenomeNode *fn,
                                                 const char *function,
                                                 GtError *err)
{
  GT_UNUSED GtNodeVisitor **vis;
  GtLuaCustomVisitor *lcv;
  GT_UNUSED GtGenomeNode **node;
  int had_err = 0;
  gt_assert(nv);
  lcv = lua_custom_visitor_cast(nv);

  node = check_genome_node(lcv->L, 1);
  vis = check_genome_visitor(lcv->L, 2);
  gt_assert(*node == (GtGenomeNode*) fn);
  gt_assert(*vis == (GtNodeVisitor*) nv);
  lua_pushvalue(lcv->L, 2);
  lua_pushstring(lcv->L, function);
  lua_gettable(lcv->L, 2);
  if (lua_isnil(lcv->L, -1)) {
    lua_pop(lcv->L, 1);
    return had_err;
  }
  lua_pushvalue(lcv->L, 2);
  gt_lua_genome_node_push(lcv->L, gt_genome_node_ref((GtGenomeNode*) fn));
  if (lua_pcall(lcv->L, 2, 0, 0)) {
    const char *error = lua_tostring(lcv->L, -1);
    gt_error_set(err, "%s", error);
    had_err = -1;
  }
  return had_err;
}

static int lua_custom_visitor_visit_feature(GtNodeVisitor *nv,
                                            GtFeatureNode *gn,
                                            GtError *err)
{
  return lua_custom_visitor_visit_node_generic(nv, (GtGenomeNode*) gn,
                                               "visit_feature", err);
}

static int lua_custom_visitor_visit_comment(GtNodeVisitor *nv,
                                            GtCommentNode *gn,
                                            GtError *err)
{
  return lua_custom_visitor_visit_node_generic(nv, (GtGenomeNode*) gn,
                                               "visit_comment", err);
}

static int lua_custom_visitor_visit_region(GtNodeVisitor *nv,
                                           GtRegionNode *gn,
                                           GtError *err)
{
  return lua_custom_visitor_visit_node_generic(nv, (GtGenomeNode*) gn,
                                               "visit_region", err);
}

static int lua_custom_visitor_visit_eof(GtNodeVisitor *nv,
                                        GtEOFNode *gn,
                                        GtError *err)
{
  return lua_custom_visitor_visit_node_generic(nv, (GtGenomeNode*) gn,
                                               "visit_eof", err);
}

static int lua_custom_visitor_visit_sequence(GtNodeVisitor *nv,
                                             GtSequenceNode *gn,
                                             GtError *err)
{
  return lua_custom_visitor_visit_node_generic(nv, (GtGenomeNode*) gn,
                                               "visit_sequence", err);
}

static int lua_custom_visitor_visit_meta(GtNodeVisitor *nv,
                                         GtMetaNode *gn,
                                         GtError *err)
{
  return lua_custom_visitor_visit_node_generic(nv, (GtGenomeNode*) gn,
                                               "visit_meta", err);
}

const GtNodeVisitorClass* lua_custom_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtLuaCustomVisitor),
                                    NULL,
                                    lua_custom_visitor_visit_comment,
                                    lua_custom_visitor_visit_feature ,
                                    lua_custom_visitor_visit_region,
                                    lua_custom_visitor_visit_sequence,
                                    lua_custom_visitor_visit_eof);
    gt_node_visitor_class_set_meta_node_func(nvc,
                                             lua_custom_visitor_visit_meta);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

static int custom_visitor_lua_index(lua_State* L)
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

static int custom_visitor_lua_newindex(lua_State* L)
{
  lua_getfenv(L, -3);
  lua_pushvalue(L, -3);
  lua_pushvalue(L, -3);
  lua_rawset(L, -3);
  return 0;
}

static int custom_visitor_lua_new(lua_State *L)
{
  GtNodeVisitor **gv;
  GtLuaCustomVisitor *lcv;
  gt_assert(L);

  gv = lua_newuserdata(L, sizeof (GtNodeVisitor*));
  gt_assert(gv);
  *gv = gt_node_visitor_create(lua_custom_visitor_class());
  gt_assert(*gv);
  lcv = lua_custom_visitor_cast(*gv);
  luaL_getmetatable(L, GENOME_VISITOR_METATABLE);
  lua_setmetatable(L, -2);

  /* set clean env for this visitor */
  lua_newtable(L);
  lua_setfenv(L, -2);
  lcv->L = L;
  return 1;
}

static int gt_node_visitor_lua_delete(lua_State *L)
{
  GtNodeVisitor **gv;
  gv = check_genome_visitor(L, 1);
  gt_node_visitor_delete(*gv);
  return 0;
}

static const struct luaL_Reg gt_node_visitor_lib_f [] = {
  { "gff3_visitor_new", gff3_visitor_lua_new },
  { "custom_visitor_new", custom_visitor_lua_new },
  { NULL, NULL }
};

int gt_lua_open_genome_visitor(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_newmetatable(L, GENOME_VISITOR_METATABLE);
  /* metatable.__index = metatable */
  /* lua_pushvalue(L, -1); */ /* duplicate the metatable */
  /* lua_setfield(L, -2, "__index"); */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, gt_node_visitor_lua_delete);
  lua_settable(L, -3);
  lua_pushstring(L, "__index");
  lua_pushcfunction(L, custom_visitor_lua_index);
  lua_settable(L, -3);
  lua_pushstring(L, "__newindex");
  lua_pushcfunction(L, custom_visitor_lua_newindex);
  lua_settable(L, -3);
  lua_pop(L, 1);
  /* register functions */
  luaL_register(L, "gt", gt_node_visitor_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}
