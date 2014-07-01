/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd.

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

#include "core/array_api.h"
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/error.h"
#include "core/gtdatapath.h"
#include "core/hashmap_api.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/str_api.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_index_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"
#include "extended/region_mapping.h"
#include "extended/spec_visitor.h"
#include "gtlua/feature_index_lua.h"
#include "gtlua/region_mapping_lua.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/gt_lua.h"
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

struct GtSpecVisitor {
  const GtNodeVisitor parent_instance;
  lua_State *L;
  GtStr *filename;
  GtHashmap *type_specs;
  int meta_ref,
      region_ref,
      comment_ref,
      sequence_ref,
      target_ref;
  GtGenomeNode *current_node;
  const char *current_aspect,
             *matcher_name;
  GtArray *graph_context;
  GtSpecResults *res;
  bool runtime_fail_hard;
};

typedef struct {
  GtFeatureNode *parent;
  GtError *err;
  GtSpecVisitor *sv;
} GtSpecVisitorNodeTraverseInfo;

const GtNodeVisitorClass* gt_spec_visitor_class();

#define spec_visitor_cast(GV)\
        gt_node_visitor_cast(gt_spec_visitor_class(), GV)

GT_UNUSED static void spec_print_context(GtSpecVisitor *sv)
{
  GtUword i;
  GtFeatureNode *fn2;
  if (!gt_array_size(sv->graph_context)) return;
  for (i = gt_array_size(sv->graph_context); i != 0; i--) {
    fn2 = *(GtFeatureNode**) gt_array_get(sv->graph_context, i-1);
    if (!fn2) {
      continue;
    }
  }
}

static const luaL_Reg spec_luasecurelibs[] = {
  {"", luaopen_base},
  {LUA_TABLIBNAME, luaopen_table},
  {LUA_STRLIBNAME, luaopen_string},
  {LUA_MATHLIBNAME, luaopen_math},
  {"gt", gt_lua_open_lib},
  {NULL, NULL}
};

static const luaL_Reg spec_luainsecurelibs[] = {
  /* These are functions affecting the system outside the Lua sandbox!
     They should only be loaded in environments where unwanted code execution
     is not a security issue! */
  {LUA_OSLIBNAME, luaopen_os},
  {LUA_IOLIBNAME, luaopen_io},
  {LUA_LOADLIBNAME, luaopen_package},
  {NULL, NULL}
};

static void spec_luaL_opencustomlibs(lua_State *L, const luaL_Reg *lib)
{
  for (; lib->func; lib++) {
    lua_pushcfunction(L, lib->func);
    lua_pushstring(L, lib->name);
    lua_call(L, 1, 0);
  }
}

static int spec_visitor_process_node(GtSpecVisitor *sv, GtFeatureNode *node,
                                     GT_UNUSED GtFeatureNode *parent,
                                     GtError *err)
{
  int had_err = 0, *ref;
  const char *type;

  if (!gt_feature_node_is_pseudo(node)) {
    type = gt_feature_node_get_type(node);
    gt_assert(type);
    if ((ref = gt_hashmap_get(sv->type_specs, type))) {
      sv->current_node = (GtGenomeNode*) node;
      lua_rawgeti(sv->L, LUA_REGISTRYINDEX, *ref);
      gt_lua_genome_node_push(sv->L, gt_genome_node_ref((GtGenomeNode*) node));
      if (lua_pcall(sv->L, 1, 0, 0)) {
        const char *error = lua_tostring(sv->L, -1);
        gt_error_set(err, "%s", error);
        had_err = -1;
      }
    }
  }

  return had_err;
}

static int spec_visitor_visit_child(GtFeatureNode* fn, void *nti, GtError *err)
{
  int had_err = 0;
  GtSpecVisitorNodeTraverseInfo* info = (GtSpecVisitorNodeTraverseInfo*) nti;

  if (gt_feature_node_has_children(fn)) {
    GtFeatureNode *oldparent = info->parent;

    had_err = spec_visitor_process_node(info->sv, fn, info->parent, err);
    if (!had_err) {
      gt_array_add(info->sv->graph_context, info->parent);
      info->parent = fn;
      had_err = gt_feature_node_traverse_direct_children(fn, info,
                                                       spec_visitor_visit_child,
                                                       err);
    }
    if (!had_err) {
      info->parent = oldparent;
    }
  } else {
    gt_array_add(info->sv->graph_context, info->parent);
    had_err = spec_visitor_process_node(info->sv, fn, info->parent, err);
  }
  if (!had_err)
    (void) gt_array_pop(info->sv->graph_context);
  return had_err;
}

static int spec_visitor_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                     GtError *err)
{
  int had_err = 0;
  GtSpecVisitor *sv;
  GtSpecVisitorNodeTraverseInfo info;
  gt_error_check(err);

  sv = spec_visitor_cast(nv);
  had_err = spec_visitor_process_node(sv, fn, NULL, err);
  if (!had_err && gt_feature_node_has_children(fn)) {
    info.sv = sv;
    info.parent = fn;
    info.err = err;
    had_err = gt_feature_node_traverse_direct_children(fn, &info,
                                                       spec_visitor_visit_child,
                                                       err);
  }
  gt_spec_results_add_cc(sv->res);
  return had_err;
}

static int spec_visitor_meta_node(GtNodeVisitor *nv, GtMetaNode *mn,
                                  GtError *err)
{
  int had_err = 0;
  GtSpecVisitor *sv;
  gt_error_check(err);

  sv = spec_visitor_cast(nv);
  if (sv->meta_ref != GT_UNDEF_INT) {
    sv->current_node = (GtGenomeNode*) mn;
    lua_rawgeti(sv->L, LUA_REGISTRYINDEX, sv->meta_ref);
    gt_lua_genome_node_push(sv->L, gt_genome_node_ref((GtGenomeNode*) mn));
    if (lua_pcall(sv->L, 1, 0, 0)) {
      const char *error = lua_tostring(sv->L, -1);
      gt_error_set(err, "%s", error);
      had_err = -1;
    }
  }
  gt_spec_results_add_cc(sv->res);
  return had_err;
}

static int spec_visitor_region_node(GtNodeVisitor *nv, GtRegionNode *rn,
                                    GtError *err)
{
  int had_err = 0;
  GtSpecVisitor *sv;
  gt_error_check(err);

  sv = spec_visitor_cast(nv);
  if (sv->region_ref != GT_UNDEF_INT) {
    sv->current_node = (GtGenomeNode*) rn;
    lua_rawgeti(sv->L, LUA_REGISTRYINDEX, sv->region_ref);
    gt_lua_genome_node_push(sv->L, gt_genome_node_ref((GtGenomeNode*) rn));
    if (lua_pcall(sv->L, 1, 0, 0)) {
      const char *error = lua_tostring(sv->L, -1);
      gt_error_set(err, "%s", error);
      had_err = -1;
    }
  }
  gt_spec_results_add_cc(sv->res);
  return had_err;
}

static int spec_visitor_comment_node(GtNodeVisitor *nv, GtCommentNode *cn,
                                     GtError *err)
{
  int had_err = 0;
  GtSpecVisitor *sv;
  gt_error_check(err);

  sv = spec_visitor_cast(nv);
  if (sv->comment_ref != GT_UNDEF_INT) {
    sv->current_node = (GtGenomeNode*) cn;
    lua_rawgeti(sv->L, LUA_REGISTRYINDEX, sv->comment_ref);
    gt_lua_genome_node_push(sv->L, gt_genome_node_ref((GtGenomeNode*) cn));
    if (lua_pcall(sv->L, 1, 0, 0)) {
      const char *error = lua_tostring(sv->L, -1);
      gt_error_set(err, "%s", error);
      had_err = -1;
    }
  }
  gt_spec_results_add_cc(sv->res);
  return had_err;
}

static int spec_visitor_sequence_node(GtNodeVisitor *nv, GtSequenceNode *sn,
                                      GtError *err)
{
  int had_err = 0;
  GtSpecVisitor *sv;
  gt_error_check(err);

  sv = spec_visitor_cast(nv);
  if (sv->sequence_ref != GT_UNDEF_INT) {
    sv->current_node = (GtGenomeNode*) sn;
    lua_rawgeti(sv->L, LUA_REGISTRYINDEX, sv->sequence_ref);
    gt_lua_genome_node_push(sv->L, gt_genome_node_ref((GtGenomeNode*) sn));
    if (lua_pcall(sv->L, 1, 0, 0)) {
      const char *error = lua_tostring(sv->L, -1);
      gt_error_set(err, "%s", error);
      had_err = -1;
    }
  }
  gt_spec_results_add_cc(sv->res);
  return had_err;
}

static void spec_visitor_free(GtNodeVisitor* nv)
{
  GtSpecVisitor *sv;
  if (!nv) return;
  sv = spec_visitor_cast(nv);
  gt_str_delete(sv->filename);
  lua_close(sv->L);
  gt_hashmap_delete(sv->type_specs);
  gt_array_delete(sv->graph_context);
}

const GtNodeVisitorClass* gt_spec_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtSpecVisitor),
                                    spec_visitor_free,
                                    spec_visitor_comment_node,
                                    spec_visitor_feature_node,
                                    spec_visitor_region_node,
                                    spec_visitor_sequence_node,
                                    NULL);
    gt_node_visitor_class_set_meta_node_func(nvc, spec_visitor_meta_node);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

/* registry keys for spec definitions */
static const char *spec_defuserdata = "gt_spec_userdata";

static int spec_register_feature_callback(lua_State *L)
{
  int ref, *refptr;
  const char *type;
  GtSpecVisitor *sv;

  /* get parameters from stack */
  ref = luaL_ref(L, LUA_REGISTRYINDEX);
  type = lua_tostring(L, -1);
  gt_log_log("registering feature specs for type %s at ref %d", type, ref);
  /* get visitor from registry */
  lua_pushlightuserdata(L, (void *) &spec_defuserdata);
  lua_gettable(L, LUA_REGISTRYINDEX);
  sv = lua_touserdata(L, -1);
  /* register Lua callback for type */
  if (!gt_hashmap_get(sv->type_specs, type)) {
    refptr = gt_malloc(sizeof (int));
    *refptr = ref;
    gt_hashmap_add(sv->type_specs, gt_cstr_dup(type), refptr);
  } else {
    luaL_where(L, 1);
    lua_pushstring(L, "duplicate definition of spec for feature type '");
    lua_pushstring(L, type);
    lua_pushstring(L, "'");
    lua_concat(L, 4);
    return lua_error(L);
  }
  return 0;
}

static int spec_register_meta_callback(lua_State *L)
{
  int ref;
  GtSpecVisitor *sv;

  /* get parameters from stack */
  ref = luaL_ref(L, LUA_REGISTRYINDEX);
  /* get visitor from registry */
  lua_pushlightuserdata(L, (void *) &spec_defuserdata);
  lua_gettable(L, LUA_REGISTRYINDEX);
  sv = lua_touserdata(L, -1);
  /* register Lua callback */
  if (sv->meta_ref == GT_UNDEF_INT) {
    sv->meta_ref = ref;
    gt_log_log("registering meta specs at ref %d", ref);
    gt_assert(sv->meta_ref != GT_UNDEF_INT);
  } else {
    luaL_where(L, 1);
    lua_pushstring(L, "duplicate definition of spec for meta nodes");
    lua_concat(L, 2);
    return lua_error(L);
  }
  return 0;
}

static int spec_register_region_callback(lua_State *L)
{
  int ref;
  GtSpecVisitor *sv;

  /* get parameters from stack */
  ref = luaL_ref(L, LUA_REGISTRYINDEX);
  /* get visitor from registry */
  lua_pushlightuserdata(L, (void *) &spec_defuserdata);
  lua_gettable(L, LUA_REGISTRYINDEX);
  sv = lua_touserdata(L, -1);
  /* register Lua callback */
  if (sv->region_ref == GT_UNDEF_INT) {
    sv->region_ref = ref;
    gt_log_log("registering region specs at ref %d", ref);
    gt_assert(sv->region_ref != GT_UNDEF_INT);
  } else {
    luaL_where(L, 1);
    lua_pushstring(L, "duplicate definition of spec for region nodes");
    lua_concat(L, 2);
    return lua_error(L);
  }
  return 0;
}

static int spec_register_comment_callback(lua_State *L)
{
  int ref;
  GtSpecVisitor *sv;

  /* get parameters from stack */
  ref = luaL_ref(L, LUA_REGISTRYINDEX);
  /* get visitor from registry */
  lua_pushlightuserdata(L, (void *) &spec_defuserdata);
  lua_gettable(L, LUA_REGISTRYINDEX);
  sv = lua_touserdata(L, -1);
  /* register Lua callback */
  if (sv->comment_ref == GT_UNDEF_INT) {
    sv->comment_ref = ref;
    gt_log_log("registering comment specs at ref %d", ref);
    gt_assert(sv->comment_ref != GT_UNDEF_INT);
  } else {
    luaL_where(L, 1);
    lua_pushstring(L, "duplicate definition of spec for comment nodes");
    lua_concat(L, 2);
    return lua_error(L);
  }
  return 0;
}

static int spec_register_sequence_callback(lua_State *L)
{
  int ref;
  GtSpecVisitor *sv;

  /* get parameters from stack */
  ref = luaL_ref(L, LUA_REGISTRYINDEX);
  /* get visitor from registry */
  lua_pushlightuserdata(L, (void *) &spec_defuserdata);
  lua_gettable(L, LUA_REGISTRYINDEX);
  sv = lua_touserdata(L, -1);
  /* register Lua callback */
  if (sv->sequence_ref == GT_UNDEF_INT) {
    sv->sequence_ref = ref;
    gt_log_log("registering sequence specs at ref %d", ref);
    gt_assert(sv->sequence_ref != GT_UNDEF_INT);
  } else {
    luaL_where(L, 1);
    lua_pushstring(L, "duplicate definition of spec for sequence nodes");
    lua_concat(L, 2);
    return lua_error(L);
  }
  return 0;
}

static int spec_it(lua_State *L)
{
  const char *name;
  GtSpecVisitor *sv;

  /* get parameters from stack */
  name = lua_tostring(L, -2);

  /* get visitor from registry */
  lua_pushlightuserdata(L, (void *) &spec_defuserdata);
  lua_gettable(L, LUA_REGISTRYINDEX);
  sv = lua_touserdata(L, -1);
  lua_pop(L, 1);

  /* handle aspect */
  sv->current_aspect = name;

  if (lua_pcall(L, 0, 0, 0)) {
    if (sv->runtime_fail_hard) {
      return lua_error(L);
    } else {
      const char *error = lua_tostring(sv->L, -1);
      gt_spec_results_add_result(sv->res, sv->current_aspect, sv->current_node,
                                 GT_SPEC_RUNTIME_ERROR, error);
    }
  }

  return 0;
}

static int spec_expect_matchdispatch(lua_State *L)
{
  GtSpecVisitor *sv;
  bool success;
  const char *msg = "";

  /* get visitor */
  lua_pushlightuserdata(L, (void *) &spec_defuserdata);
  lua_gettable(L, LUA_REGISTRYINDEX);
  sv = lua_touserdata(L, -1);

  /* get correct matcher function from speclib */
  lua_getglobal(L, "matchers");
  lua_pushstring(L, sv->matcher_name);
  lua_gettable(L, -2);
  lua_insert(L, 1);

  /* get target */
  lua_rawgeti(L, LUA_REGISTRYINDEX, sv->target_ref);
  lua_insert(L, 2);

  /* call matcher */
  lua_pop(L, 2);
  if (!lua_isfunction(L, 1)) {
    luaL_where(L, 1);
    lua_pushstring(L, "matcher '");
    lua_pushstring(L, sv->matcher_name);
    lua_pushstring(L, "' is not a function");
    lua_concat(L, 4);
    return lua_error(L);
  }
  lua_pcall(L, lua_gettop(L)-1, 2, 0);

  /* handle results */
  if (!lua_isboolean(L, 1)) {
    lua_pushstring(L, "matcher '");
    lua_pushstring(L, sv->matcher_name);
    lua_pushstring(L, "' must return boolean");
    lua_concat(L, 3);
    return lua_error(L);
  }
  success = lua_toboolean(L, 1);
  if (!success) {
    gt_assert(lua_isstring(L, 2));
    msg = lua_tostring(L, 2);
    luaL_where(L, 1);
    lua_pushstring(L, msg);
    lua_concat(L, 2);
    msg = lua_tostring(L, -1);
  }
  gt_spec_results_add_result(sv->res, sv->current_aspect, sv->current_node,
                             success ? GT_SPEC_SUCCESS : GT_SPEC_FAILURE, msg);

  return 0;
}

static int spec_expect_index(lua_State *L)
{
  const char *matcher_name;
  GtSpecVisitor *sv;

  /* get parameters from stack */
  matcher_name = lua_tostring(L, 2);
  lua_pushlightuserdata(L, (void *) &spec_defuserdata);
  lua_gettable(L, LUA_REGISTRYINDEX);
  sv = lua_touserdata(L, -1);

  /* we can't do closures as nicely as in Lua -- keep this around for later */
  sv->matcher_name = matcher_name;

  /* return dispatcher function */
  lua_pushcfunction(L, spec_expect_matchdispatch);

  return 1;
}

static int spec_expect(lua_State *L)
{
  int ref;
  GtSpecVisitor *sv;

  /* accept only one parameter */
  if (lua_gettop(L) > 1) {
    luaL_where(L, 1);
    lua_pushstring(L, "'expect' takes only one parameter");
    lua_concat(L, 2);
    return lua_error(L);
  }

  /* get parameters from stack */
  ref = luaL_ref(L, LUA_REGISTRYINDEX);

  /* get visitor */
  lua_pushlightuserdata(L, (void *) &spec_defuserdata);
  lua_gettable(L, LUA_REGISTRYINDEX);
  sv = lua_touserdata(L, -1);

  /* we can't do closures as nicely as in Lua -- keep this around for later */
  sv->target_ref = ref;

  /* return new object with __index metamethod set */
  lua_newtable(L);
  luaL_getmetatable(L, "GenomeTools.specexpect");
  lua_setmetatable(L, -2);

  return 1;
}

static int spec_feature_node_lua_appears_as_child_of_type(lua_State *L)
{
  GtGenomeNode **gn;
  GtFeatureNode *fn, *fn2;
  bool found = false;
  GtSpecVisitor *sv;
  const char *type;
  GtUword i;
  gn = check_genome_node(L, 1);
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  type = gt_symbol(luaL_checkstring(L, 2));

  lua_pushlightuserdata(L, (void *) &spec_defuserdata);
  lua_gettable(L, LUA_REGISTRYINDEX);
  sv = lua_touserdata(L, -1);

  gt_assert(sv && sv->graph_context);
  if (gt_array_size(sv->graph_context) > 0) {
    for (i = gt_array_size(sv->graph_context) - 1; i != 0; i--) {
      fn2 = *(GtFeatureNode**) gt_array_get(sv->graph_context, i);
      if (!fn2) continue;
      found = (gt_feature_node_get_type(fn2) == type);
    }
  }
  lua_pushboolean(L, found);
  return 1;
}

static int spec_init_lua_env(GtSpecVisitor *sv, const char *progname)
{
  int had_err = 0;
  GtStr *prog, *speclib;
  gt_assert(sv);

  /* add some more convenience methods not in the public API */
  luaL_getmetatable(sv->L, GENOME_NODE_METATABLE);
  lua_pushstring(sv->L, "appears_as_child_of_type");
  lua_pushcfunction(sv->L, spec_feature_node_lua_appears_as_child_of_type);
  lua_rawset(sv->L, -3);

  /* setup DSL: node-type-specific 'describe' environment */
  lua_newtable(sv->L);
  lua_pushstring(sv->L, "feature");
  lua_pushcfunction(sv->L, spec_register_feature_callback);
  lua_rawset(sv->L, -3);
  lua_pushstring(sv->L, "region");
  lua_pushcfunction(sv->L, spec_register_region_callback);
  lua_rawset(sv->L, -3);
  lua_pushstring(sv->L, "meta");
  lua_pushcfunction(sv->L, spec_register_meta_callback);
  lua_rawset(sv->L, -3);
  lua_pushstring(sv->L, "comment");
  lua_pushcfunction(sv->L, spec_register_comment_callback);
  lua_rawset(sv->L, -3);
  lua_pushstring(sv->L, "sequence");
  lua_pushcfunction(sv->L, spec_register_sequence_callback);
  lua_rawset(sv->L, -3);
  lua_setglobal(sv->L, "describe");

  /* setup DSL: aspect definition 'it' */
  lua_pushcfunction(sv->L, spec_it);
  lua_setglobal(sv->L, "it");

  /* setup DSL: 'expect'ations */
  lua_pushcfunction(sv->L, spec_expect);
  lua_setglobal(sv->L, "expect");

  /* setup DSL: dispatcher metatable for Lua-based expectation matchers */
  luaL_newmetatable(sv->L, "GenomeTools.specexpect");
  lua_pushstring(sv->L, "__index");
  lua_pushcfunction(sv->L, spec_expect_index);
  lua_settable(sv->L, -3);

  /* load matcher funcs written in Lua (for extensibility) */
  prog = gt_str_new();
  gt_str_append_cstr_nt(prog, progname,
                        gt_cstr_length_up_to_char(progname, ' '));
  speclib = gt_get_gtdata_path(gt_str_get(prog), NULL);
  gt_str_delete(prog);
  gt_str_append_cstr(speclib, "/spec/speclib.lua");
  had_err = (luaL_loadfile(sv->L, gt_str_get(speclib))
               || lua_pcall(sv->L, 0, 0, 0));
  gt_assert(!had_err);
  gt_str_delete(speclib);

  /* store this visitor for later use */
  lua_pushlightuserdata(sv->L, (void*) &spec_defuserdata);
  lua_pushlightuserdata(sv->L, (void*) sv);
  lua_settable(sv->L, LUA_REGISTRYINDEX);

  /* no meta and region specs at the beginning */
  sv->meta_ref = sv->region_ref =
    sv->sequence_ref = sv->comment_ref = GT_UNDEF_INT;

  return had_err;
}

void gt_spec_visitor_report_runtime_errors(GtSpecVisitor *sv)
{
  gt_assert(sv);
  sv->runtime_fail_hard = false;
}

void gt_spec_visitor_fail_on_runtime_error(GtSpecVisitor *sv)
{
  gt_assert(sv);
  sv->runtime_fail_hard = true;
}

void gt_spec_visitor_add_feature_index(GtSpecVisitor *sv, GtFeatureIndex *fi)
{
  gt_assert(sv && sv->L && fi);
  gt_lua_feature_index_push(sv->L, fi);
  lua_setglobal(sv->L, "feature_index");
}

void gt_spec_visitor_add_region_mapping(GtSpecVisitor *sv, GtRegionMapping *rm)
{
  gt_assert(sv && sv->L && rm);
  gt_lua_region_mapping_push(sv->L, rm);
  lua_setglobal(sv->L, "region_mapping");
}

GtNodeVisitor* gt_spec_visitor_new(const char *specfile, GtSpecResults *res,
                                   GtError *err)
{
  GtNodeVisitor *nv;
  GtSpecVisitor *sv;
  gt_assert(specfile);
  nv = gt_node_visitor_create(gt_spec_visitor_class());
  sv = spec_visitor_cast(nv);
  sv->L = luaL_newstate();
  if (!sv->L) {
    gt_error_set(err, "cannot create new Lua state");
    gt_node_visitor_delete(nv);
    return NULL;
  }
  spec_luaL_opencustomlibs(sv->L, spec_luasecurelibs);
  spec_luaL_opencustomlibs(sv->L, spec_luainsecurelibs);
  sv->filename = gt_str_new_cstr(specfile);
  sv->res = res;
  sv->type_specs = gt_hashmap_new(GT_HASH_STRING, gt_free_func, gt_free_func);
  sv->graph_context = gt_array_new(sizeof (GtFeatureNode*));
  sv->runtime_fail_hard = false;

  (void) spec_init_lua_env(sv, gt_error_get_progname(err));
  if (luaL_loadfile(sv->L, specfile) || lua_pcall(sv->L, 0, 0, 0)) {
    gt_error_set(err, "%s", lua_tostring(sv->L, -1));
    gt_node_visitor_delete(nv);
    return NULL;
  }
  return nv;
}
