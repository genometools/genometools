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

#include <stdio.h>
#include "core/array.h"
#include "core/cstr_api.h"
#include "core/file.h"
#include "core/hashmap_api.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/str.h"
#include "core/str_array_api.h"
#include "core/symbol_api.h"
#include "core/unused_api.h"
#include "extended/comment_node_api.h"
#include "extended/feature_node_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_visitor.h"
#include "extended/luahelper.h"
#include "extended/spec_results.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/genome_visitor_lua.h"
#include "gtlua/gt_lua.h"
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

typedef struct {
  GtStrArray *failure_messages;
  GtStrArray *runtime_error_messages;
} GtSpecAspectNodeResult;

typedef struct {
  GtUword *fail, *err;
} GtSpecAspectCountInfo;

typedef struct {
  GtHashmap *aspect_node_results;
  GtStr *name;
  GtGenomeNode *last_node;  /* not thread-safe */
  GtUword nof_nodes;
} GtSpecAspect;

typedef struct {
  GtFile *outfile;
  bool details,
       colored;
  GtUword node_i;
  lua_State *L;
  GtSpecResults *sr;
} GtSpecResultsReportInfo;

struct GtSpecResults {
  GtHashmap *feature_aspects,
            *meta_aspects,
            *region_aspects,
            *comment_aspects,
            *sequence_aspects;
  bool seen_feature,
       seen_meta,
       seen_region,
       seen_comment,
       seen_sequence,
       has_failures,
       has_runtime_errors;
  GtUword checked_types,
          checked_aspects,
          checked_ccs,
          checked_nodes;
  GtStrArray *warnings;
};

static const luaL_Reg spec_results_luasecurelibs[] = {
  {"", luaopen_base},
  {LUA_TABLIBNAME, luaopen_table},
  {LUA_STRLIBNAME, luaopen_string},
  {LUA_MATHLIBNAME, luaopen_math},
  {LUA_IOLIBNAME, luaopen_io},
  {LUA_OSLIBNAME, luaopen_os},
  {LUA_LOADLIBNAME, luaopen_package},
  {"gt", gt_lua_open_lib},
  {NULL, NULL}
};

static void spec_results_luaL_opencustomlibs(lua_State *L, const luaL_Reg *lib)
{
  for (; lib->func; lib++) {
    lua_pushcfunction(L, lib->func);
    lua_pushstring(L, lib->name);
    lua_call(L, 1, 0);
  }
}

static GtSpecAspectNodeResult* gt_spec_aspect_node_result_new()
{
  GtSpecAspectNodeResult *sanr;

  sanr = gt_calloc(1, sizeof (*sanr));
  sanr->failure_messages = gt_str_array_new();
  sanr->runtime_error_messages = gt_str_array_new();
  return sanr;
}

static void gt_spec_aspect_node_result_delete(GtSpecAspectNodeResult *sanr)
{
  if (!sanr) return;
  gt_str_array_delete(sanr->failure_messages);
  gt_str_array_delete(sanr->runtime_error_messages);
  gt_free(sanr);
}

static GtSpecAspect* gt_spec_aspect_new(const char *name)
{
  GtSpecAspect *sa;
  gt_assert(name);

  sa = gt_calloc(1, sizeof (*sa));
  sa->name = gt_str_new_cstr(name);
  sa->last_node = NULL;
  sa->aspect_node_results= gt_hashmap_new(GT_HASH_DIRECT,
                                    (GtFree) gt_genome_node_delete,
                                    (GtFree) gt_spec_aspect_node_result_delete);

  return sa;
}

static void gt_spec_aspect_delete(GtSpecAspect *sa)
{
  if (!sa) return;
  gt_str_delete(sa->name);
  gt_hashmap_delete(sa->aspect_node_results);
  gt_free(sa);
}

GtSpecResults *gt_spec_results_new(void)
{
  GtSpecResults *spec_results;
  spec_results = gt_calloc(1, sizeof (GtSpecResults));
  spec_results->feature_aspects = gt_hashmap_new(GT_HASH_STRING, NULL,
                                                 (GtFree) gt_hashmap_delete);
  spec_results->meta_aspects = gt_hashmap_new(GT_HASH_STRING, NULL,
                                              (GtFree) gt_spec_aspect_delete);
  spec_results->sequence_aspects = gt_hashmap_new(GT_HASH_STRING, NULL,
                                               (GtFree) gt_spec_aspect_delete);
  spec_results->region_aspects = gt_hashmap_new(GT_HASH_STRING, NULL,
                                               (GtFree) gt_spec_aspect_delete);
  spec_results->comment_aspects = gt_hashmap_new(GT_HASH_STRING, NULL,
                                               (GtFree) gt_spec_aspect_delete);
  spec_results->warnings = gt_str_array_new();
  spec_results->has_runtime_errors = spec_results->has_failures = false;
  return spec_results;
}

void gt_spec_results_add_cc(GtSpecResults *sr)
{
  gt_assert(sr);
  sr->checked_ccs++;
}

void gt_spec_results_record_warning(GtSpecResults *sr, const char *w)
{
  gt_assert(sr && w);
  gt_str_array_add_cstr(sr->warnings, w);
}

void gt_spec_results_add_result(GtSpecResults *sr,
                                const char *aspect,
                                GtGenomeNode *node,
                                GtSpecResultStatus status,
                                const char *error_string)
{
  GtSpecAspect *sa = NULL;
  GtSpecAspectNodeResult *sanr = NULL;
  gt_assert(aspect && node && error_string);

  if (gt_feature_node_try_cast(node)) {
    GtHashmap *per_type_aspects;
    const char *type = gt_feature_node_get_type((GtFeatureNode*) node);
    if (!(per_type_aspects = gt_hashmap_get(sr->feature_aspects, type))) {
      per_type_aspects = gt_hashmap_new(GT_HASH_STRING, NULL,
                                        (GtFree) gt_spec_aspect_delete);
      gt_hashmap_add(sr->feature_aspects, (void*) gt_symbol(type),
                     per_type_aspects);
      sr->checked_types++;
      sr->seen_feature = true;
    }
    gt_assert(per_type_aspects);
    if (!(sa = gt_hashmap_get(per_type_aspects, aspect))) {
      sa = gt_spec_aspect_new(aspect);
      gt_hashmap_add(per_type_aspects, (void*) gt_symbol(aspect), sa);
      sr->checked_aspects++;
    }
  } else if (gt_meta_node_try_cast(node)) {
    if (!(sa = gt_hashmap_get(sr->meta_aspects, aspect))) {
      sa = gt_spec_aspect_new(aspect);
      gt_hashmap_add(sr->meta_aspects, (void*) gt_symbol(aspect), sa);
      sr->checked_aspects++;
      sr->seen_meta = true;
    }
  } else if (gt_region_node_try_cast(node)) {
    if (!(sa = gt_hashmap_get(sr->region_aspects, aspect))) {
      sa = gt_spec_aspect_new(aspect);
      gt_hashmap_add(sr->region_aspects, (void*) gt_symbol(aspect), sa);
      sr->checked_aspects++;
      sr->seen_region = true;
    }
  } else if (gt_comment_node_try_cast(node)) {
    if (!(sa = gt_hashmap_get(sr->comment_aspects, aspect))) {
      sa = gt_spec_aspect_new(aspect);
      gt_hashmap_add(sr->comment_aspects, (void*) gt_symbol(aspect), sa);
      sr->checked_aspects++;
      sr->seen_comment = true;
    }
  } else if (gt_sequence_node_try_cast(node)) {
    if (!(sa = gt_hashmap_get(sr->sequence_aspects, aspect))) {
      sa = gt_spec_aspect_new(aspect);
      gt_hashmap_add(sr->sequence_aspects, (void*) gt_symbol(aspect), sa);
      sr->checked_aspects++;
      sr->seen_sequence = true;
    }
  }
  gt_assert(sa);
  /* count number of nodes seen for this aspect */
  if (sa->last_node != node)
    sa->nof_nodes++;
  switch (status) {
    case GT_SPEC_SUCCESS:
      /* pass */
      break;
    case GT_SPEC_FAILURE:
      if (!(sanr = gt_hashmap_get(sa->aspect_node_results, node))) {
        sanr = gt_spec_aspect_node_result_new();
        node = gt_genome_node_ref(node);
        gt_hashmap_add(sa->aspect_node_results, node, sanr);
      }
      gt_assert(sanr);
      gt_str_array_add_cstr(sanr->failure_messages, error_string);
      if (!sr->has_failures)
        sr->has_failures = true;
      break;
    case GT_SPEC_RUNTIME_ERROR:
      if (!(sanr = gt_hashmap_get(sa->aspect_node_results, node))) {
        sanr = gt_spec_aspect_node_result_new();
        node = gt_genome_node_ref(node);
        gt_hashmap_add(sa->aspect_node_results, node, sanr);
      }
      gt_assert(sanr);
      gt_str_array_add_cstr(sanr->runtime_error_messages, error_string);
      if (!sr->has_runtime_errors)
        sr->has_runtime_errors = true;
      break;
  }
  sa->last_node = node;
}

static int gt_spec_aspect_make_node_model(void *key, void *value, void *data,
                                          GT_UNUSED GtError *err)
{
  GtGenomeNode *gn = (GtGenomeNode*) key;
  GtSpecAspectNodeResult *sanr = (GtSpecAspectNodeResult*) value;
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  const char *tmp = NULL;
  bool has_id = false;
  GtUword i = 0;

  if (gt_feature_node_try_cast(gn)) {
    if ((tmp = gt_feature_node_get_attribute((GtFeatureNode*) gn, "ID"))) {
      has_id = true;
    }
  }
  if (gt_str_array_size(sanr->failure_messages) > 0
        || gt_str_array_size(sanr->runtime_error_messages) > 0) {
    lua_pushnumber(info->L, info->node_i++);
    lua_newtable(info->L);
    if (has_id) {
      lua_pushstring(info->L, "ID");
      lua_pushstring(info->L, tmp);
      lua_rawset(info->L, -3);
    }
    lua_pushstring(info->L, "filename");
    lua_pushstring(info->L, gt_genome_node_get_filename(gn));
    lua_rawset(info->L, -3);
    lua_pushstring(info->L, "linenumber");
    lua_pushnumber(info->L, gt_genome_node_get_line_number(gn));
    lua_rawset(info->L, -3);
    lua_pushstring(info->L, "node");
    gt_lua_genome_node_push(info->L, gt_genome_node_ref(gn));
    lua_rawset(info->L, -3);
    lua_pushstring(info->L, "failure_messages");
    lua_newtable(info->L);
    for (i = 0; i < gt_str_array_size(sanr->failure_messages); i++) {
      lua_pushnumber(info->L, i+1);
      lua_pushstring(info->L, gt_str_array_get(sanr->failure_messages, i));
      lua_rawset(info->L, -3);
    }
    lua_rawset(info->L, -3);
    lua_pushstring(info->L, "runtime_error_messages");
    lua_newtable(info->L);
    for (i = 0; i < gt_str_array_size(sanr->runtime_error_messages); i++) {
      lua_pushnumber(info->L, i+1);
      lua_pushstring(info->L, gt_str_array_get(sanr->runtime_error_messages,
                                               i));
      lua_rawset(info->L, -3);
    }
    lua_rawset(info->L, -3);
    lua_rawset(info->L, -3);
  }
  return 0;
}

static int gt_spec_aspect_count_stats(GT_UNUSED void *key, void *value,
                                      void *data, GT_UNUSED GtError *err)
{
  GtSpecAspectNodeResult *sanr = (GtSpecAspectNodeResult*) value;
  GtSpecAspectCountInfo *info = (GtSpecAspectCountInfo*) data;

  if (gt_str_array_size(sanr->runtime_error_messages) > 0) {
    (*info->err)++;
  } else if (gt_str_array_size(sanr->failure_messages) > 0) {
    (*info->fail)++;
  }
  return 0;
}

static int gt_spec_results_make_aspect_model(GT_UNUSED void *key, void *value,
                                             void *data, GT_UNUSED GtError *err)
{
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  GtSpecAspect* sa =  (GtSpecAspect*) value;
  GtSpecAspectCountInfo count_info;
  GT_UNUSED int rval;
  GtUword successes = sa->nof_nodes,
          failures = 0,
          runtime_errors = 0;
  gt_assert(sa && info);

  /* update global number of checked node/aspect pairs */
  info->sr->checked_nodes += sa->nof_nodes;

  /* calculate number of failures/errors */
  count_info.fail = &failures;
  count_info.err = &runtime_errors;
  rval = gt_hashmap_foreach(sa->aspect_node_results, gt_spec_aspect_count_stats,
                            &count_info, NULL);
  gt_assert(!rval); /* gt_spec_aspect_count_stats() is sane */

  /* calculate number of successes */
  gt_assert(failures <= sa->nof_nodes && runtime_errors <= sa->nof_nodes);
  successes = sa->nof_nodes - failures - runtime_errors;
  gt_assert(failures + runtime_errors + successes == sa->nof_nodes);

  /* export results as data model */
  info->node_i = 1;
  lua_pushstring(info->L, gt_str_get(sa->name));
  lua_newtable(info->L);
  lua_pushstring(info->L, "successes");
  lua_pushnumber(info->L, successes);
  lua_rawset(info->L, -3);
  lua_pushstring(info->L, "failures");
  lua_pushnumber(info->L, failures);
  lua_rawset(info->L, -3);
  lua_pushstring(info->L, "runtime_errors");
  lua_pushnumber(info->L, runtime_errors);
  lua_rawset(info->L, -3);
  lua_pushstring(info->L, "nodes");
  lua_newtable(info->L);
  rval = gt_hashmap_foreach(sa->aspect_node_results,
                              gt_spec_aspect_make_node_model, info, NULL);
  gt_assert(rval == 0);
  lua_rawset(info->L, -3);
  lua_rawset(info->L, -3);

  return 0;
}

static int gt_spec_results_make_feature_model(void *key, void *value,
                                              void *data,
                                                         GT_UNUSED GtError *err)
{
  const char *type = (const char*) key;
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  lua_pushstring(info->L, type);
  lua_newtable(info->L);
  gt_hashmap_foreach((GtHashmap*) value, gt_spec_results_make_aspect_model,
                     info, NULL);
  lua_rawset(info->L, -3);
  return 0;
}

/* registry keys for outfile pointer */
static const char *spec_resuserdata = "gt_spec_results_userdata";

static int gt_spec_results_lua_print(lua_State* L) {
    GtFile *outfile;
    int i = 0;
    gt_assert(L);

    lua_pushlightuserdata(L, (void *) &spec_resuserdata);
    lua_gettable(L, LUA_REGISTRYINDEX);
    outfile = lua_touserdata(L, -1);

    for (i = 1; i <= lua_gettop(L); i++) {
      if (lua_isstring(L, i)) {
        gt_file_xprintf(outfile, "%s", lua_tostring(L, i));
      }
    }

    return 0;
}

bool gt_spec_results_has_runtime_errors(GtSpecResults *sr)
{
  gt_assert(sr);
  return sr->has_runtime_errors;
}

bool gt_spec_results_has_failures(GtSpecResults *sr)
{
  gt_assert(sr);
  return sr->has_failures;
}

int gt_spec_results_render_template(GtSpecResults *sr, const char *template,
                                    GtFile *outfile, const char *specfile,
                                    bool details, bool colored,
                                    const char *runtime_str, GtError *err)
{
  lua_State *L;
  GtNodeVisitor *vis = NULL;
  GtSpecResultsReportInfo info;
  int had_err = 0;
  gt_assert(sr && specfile && template && err);

  L = luaL_newstate();
  if (!L) {
    gt_error_set(err, "cannot create new Lua state");
    return -1;
  }
  spec_results_luaL_opencustomlibs(L, spec_results_luasecurelibs);

  /* fill info */
  info.outfile = outfile;
  info.details = details;
  info.colored = colored;
  info.L = L;
  info.sr = sr;

  /* make model for nodes/aspects */
  lua_newtable(L);
  gt_hashmap_foreach_in_key_order(sr->feature_aspects,
                                  gt_spec_results_make_feature_model,
                                  &info, NULL);
  lua_setglobal(L, "features");
  lua_newtable(L);
  gt_hashmap_foreach_in_key_order(sr->region_aspects,
                                  gt_spec_results_make_aspect_model,
                                  &info, NULL);
  lua_setglobal(L, "regions");
  lua_newtable(L);
  gt_hashmap_foreach_in_key_order(sr->meta_aspects,
                                  gt_spec_results_make_aspect_model,
                                  &info, NULL);
  lua_setglobal(L, "metas");
  lua_newtable(L);
  gt_hashmap_foreach_in_key_order(sr->sequence_aspects,
                                  gt_spec_results_make_aspect_model,
                                  &info, NULL);
  lua_setglobal(L, "sequences");
  lua_newtable(L);
  gt_hashmap_foreach_in_key_order(sr->comment_aspects,
                                  gt_spec_results_make_aspect_model,
                                  &info, NULL);
  lua_setglobal(L, "comments");

  /* add warnings */
  lua_newtable(L);
  if (gt_str_array_size(sr->warnings) > 0) {
    GtUword i;
    for (i = 0; i < gt_str_array_size(sr->warnings); i++) {
      lua_pushnumber(L, i+1);
      lua_pushstring(L, gt_str_array_get(sr->warnings, i));
      lua_rawset(L, -3);
    }
  }
  lua_setglobal(L, "warnings");

  /* dome global settings */
  lua_newtable(L);
  lua_pushstring(L, "runtime");
  lua_pushstring(L, runtime_str);
  lua_rawset(L, -3);
  lua_pushstring(L, "spec_filename");
  lua_pushstring(L, specfile);
  lua_rawset(L, -3);
  lua_pushstring(L, "template_filename");
  lua_pushstring(L, template);
  lua_rawset(L, -3);
  lua_pushstring(L, "show_details");
  lua_pushboolean(L, details);
  lua_rawset(L, -3);
  lua_pushstring(L, "colored_output");
  lua_pushboolean(L, colored);
  lua_rawset(L, -3);
  lua_pushstring(L, "coloured_output");
  lua_pushboolean(L, colored);
  lua_rawset(L, -3);
  lua_pushstring(L, "checked_ccs");
  lua_pushnumber(L, sr->checked_ccs);
  lua_rawset(L, -3);
  lua_pushstring(L, "checked_aspects");
  lua_pushnumber(L, sr->checked_aspects);
  lua_rawset(L, -3);
  lua_pushstring(L, "checked_nodes");
  lua_pushnumber(L, sr->checked_nodes);
  lua_rawset(L, -3);
  lua_pushstring(L, "checked_types");
  lua_pushnumber(L, sr->checked_types);
  lua_rawset(L, -3);
  lua_setglobal(L, "global");

  /* add own print function to support -o etc. */
  lua_pushcfunction(L, gt_spec_results_lua_print);
  lua_setglobal(L, "template_print");

  /* gt.script_dir is useful in templates too */
  gt_lua_set_script_dir(L, template);

  /* also push a GtGFF3Visitor writing to outfile*/
  vis = gt_gff3_visitor_new(outfile);
  gt_gff3_visitor_retain_id_attributes((GtGFF3Visitor*) vis);
  gt_gff3_visitor_allow_nonunique_ids((GtGFF3Visitor*) vis);
  gt_lua_genome_visitor_push(L, vis);
  lua_setglobal(L, "gff3_out_visitor");

  lua_pushlightuserdata(L, (void*) &spec_resuserdata);
  lua_pushlightuserdata(L, (void*) outfile);
  lua_settable(L, LUA_REGISTRYINDEX);

  /* render template */
  had_err = (luaL_loadfile(L, template)
               || lua_pcall(L, 0, 0, 0));
  if (had_err) {
    const char *error = lua_tostring(L, -1);
    gt_error_set(err, "%s", error);
    had_err = -1;
  }

  lua_close(L);
  return had_err;
}

void gt_spec_results_delete(GtSpecResults *spec_results)
{
  if (!spec_results) return;
  gt_hashmap_delete(spec_results->feature_aspects);
  gt_hashmap_delete(spec_results->meta_aspects);
  gt_hashmap_delete(spec_results->comment_aspects);
  gt_hashmap_delete(spec_results->sequence_aspects);
  gt_hashmap_delete(spec_results->region_aspects);
  gt_str_array_delete(spec_results->warnings);
  gt_free(spec_results);
}
