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
#include "extended/luahelper.h"
#include "extended/spec_results.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/gt_lua.h"
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

#define GT_SPEC_ANSI_COLOR_RED     "\x1b[31m"
#define GT_SPEC_ANSI_COLOR_GREEN   "\x1b[32m"
#define GT_SPEC_ANSI_COLOR_YELLOW  "\x1b[33m"
#define GT_SPEC_ANSI_COLOR_BLUE    "\x1b[34m"
#define GT_SPEC_ANSI_COLOR_MAGENTA "\x1b[35m"
#define GT_SPEC_ANSI_COLOR_CYAN    "\x1b[36m"
#define GT_SPEC_ANSI_COLOR_RESET   "\x1b[0m"

typedef struct {
  GtStrArray *failure_messages;
  GtStrArray *runtime_error_messages;
} GtSpecAspectNodeResult;

typedef struct {
  GtUword *succ, *fail, *err;
} GtSpecAspectCountInfo;

typedef struct {
  GtHashmap *aspect_node_results;
  GtStr *name;
  GtGenomeNode *last_node;
  GtUword node_successes,
          node_failures,
          node_runtime_errors,
          successes,
          failures,
          runtime_errors;
} GtSpecAspect;

typedef struct {
  GtFile *outfile;
  bool details,
       colored;
  GtUword node_i;
  lua_State *L;
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
       seen_sequence;
  GtUword checked_types,
          checked_aspects,
          checked_ccs,
          checked_nodes;
  GtStrArray *warnings;
  bool record_per_node; /* not thread-safe, requires sequential processing
                           for one node! */
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

static int gt_spec_aspect_report_node(void *key, void *value, void *data,
                                      GT_UNUSED GtError *err)
{
  GtGenomeNode *gn = (GtGenomeNode*) key;
  GtSpecAspectNodeResult *sanr = (GtSpecAspectNodeResult*) value;
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  char tmpbuf[BUFSIZ];
  bool has_id = false;
  GtUword i = 0;

  if (gt_feature_node_try_cast(gn)) {
    const char *tmp;
    if ((tmp = gt_feature_node_get_attribute((GtFeatureNode*) gn, "ID"))) {
      has_id = true;
      (void) snprintf(tmpbuf, BUFSIZ, "%s, ", tmp);
    }
  }
  if (gt_str_array_size(sanr->failure_messages) > 0
        || gt_str_array_size(sanr->runtime_error_messages) > 0) {
    gt_file_xprintf(info->outfile,
                    "%s      offending node #" GT_WU " "
                    "(%sfrom %s, line %u):%s\n",
                    info->colored ? GT_SPEC_ANSI_COLOR_RED : "",
                    info->node_i++,
                    (has_id ? tmpbuf : ""),
                    gt_genome_node_get_filename(gn),
                    gt_genome_node_get_line_number(gn),
                    info->colored ? GT_SPEC_ANSI_COLOR_RESET : "");
  }
  for (i = 0; i < gt_str_array_size(sanr->failure_messages); i++) {
    gt_file_xprintf(info->outfile,
                    "%s         %s%s\n",
                    info->colored ? GT_SPEC_ANSI_COLOR_RED : "",
                    gt_str_array_get(sanr->failure_messages, i),
                    info->colored ? GT_SPEC_ANSI_COLOR_RESET : "");
  }
  for (i = 0; i < gt_str_array_size(sanr->runtime_error_messages); i++) {
    gt_file_xprintf(info->outfile,
                    "%s         %s%s\n",
                    info->colored ? GT_SPEC_ANSI_COLOR_MAGENTA : "",
                    gt_str_array_get(sanr->runtime_error_messages, i),
                    info->colored ? GT_SPEC_ANSI_COLOR_RESET : "");
  }

  return 0;
}

static void gt_spec_aspect_report(GtSpecAspect *sa, GtFile *outfile,
                                  GtSpecResultsReportInfo *info)
{
  GT_UNUSED int rval;
  gt_assert(sa && info);

  /* output stats */
  info->node_i = 1;
  GtStr *outbuf = gt_str_new_cstr("  - ");
  gt_str_append_str(outbuf, sa->name);
  gt_str_append_cstr(outbuf, " (");
  if (sa->successes > 0) {
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_GREEN);
    gt_str_append_ulong(outbuf, sa->successes);
    gt_str_append_cstr(outbuf, " success");
    if (sa->successes > 1)
      gt_str_append_cstr(outbuf, "es");
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_RESET);
  }
  if (sa->failures > 0) {
    if (sa->successes > 0)
      gt_str_append_cstr(outbuf, ", ");
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_RED);
    gt_str_append_ulong(outbuf, sa->failures);
    gt_str_append_cstr(outbuf, " failure");
    if (sa->failures > 1)
      gt_str_append_cstr(outbuf, "s");
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_RESET);
  }
  if (sa->runtime_errors > 0) {
    if (sa->successes > 0 || sa->failures > 0)
      gt_str_append_cstr(outbuf, ", ");
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_MAGENTA);
    gt_str_append_ulong(outbuf, sa->runtime_errors);
    gt_str_append_cstr(outbuf, " runtime error");
    if (sa->runtime_errors > 1)
      gt_str_append_cstr(outbuf, "s");
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_RESET);
  }
  gt_str_append_cstr(outbuf, ")\n");
  gt_file_xprintf(outfile, "%s", gt_str_get(outbuf));
  gt_str_delete(outbuf);
  if (info->details) {
    rval = gt_hashmap_foreach(sa->aspect_node_results,
                              gt_spec_aspect_report_node, info, NULL);
    gt_assert(rval == 0);
  }
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
  spec_results->record_per_node = false;
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

void gt_spec_results_record_per_node(GtSpecResults *sr)
{
  gt_assert(sr);
  sr->record_per_node = true;
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
  switch (status) {
    case GT_SPEC_SUCCESS:
      if (!sr->record_per_node || sa->last_node != node) {
        sa->successes++;
      }
      break;
    case GT_SPEC_FAILURE:
      if (!sr->record_per_node || sa->last_node != node) {
        sa->failures++;
      }
      if (!(sanr = gt_hashmap_get(sa->aspect_node_results, node))) {
        sanr = gt_spec_aspect_node_result_new();
        node = gt_genome_node_ref(node);
        gt_hashmap_add(sa->aspect_node_results, node, sanr);
      }
      gt_assert(sanr);
      gt_str_array_add_cstr(sanr->failure_messages, error_string);
      break;
    case GT_SPEC_RUNTIME_ERROR:
      if (!sr->record_per_node || sa->last_node != node) {
        sa->runtime_errors++;
      }
      if (!(sanr = gt_hashmap_get(sa->aspect_node_results, node))) {
        sanr = gt_spec_aspect_node_result_new();
        node = gt_genome_node_ref(node);
        gt_hashmap_add(sa->aspect_node_results, node, sanr);
      }
      gt_assert(sanr);
      gt_str_array_add_cstr(sanr->runtime_error_messages, error_string);
      break;
  }
  sa->last_node = node;
  sr->checked_nodes++;
}

static int gt_spec_results_report_single(GT_UNUSED void *key, void *value,
                                         void *data, GT_UNUSED GtError *err)
{
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  gt_spec_aspect_report((GtSpecAspect*) value, info->outfile, info);
  return 0;
}

static int gt_spec_results_report_features(void *key, void *value, void *data,
                                           GT_UNUSED GtError *err)
{
  const char *type = (const char*) key;
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  gt_file_xprintf(info->outfile, "a %s%s%s feature\n",
                           info->colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                           type,
                           info->colored ? GT_SPEC_ANSI_COLOR_RESET : "");
  gt_hashmap_foreach((GtHashmap*) value, gt_spec_results_report_single,
                     info, NULL);
  return 0;
}

static int gt_spec_aspect_make_node_model(void *key, void *value, void *data,
                                          GT_UNUSED GtError *err)
{
  GtGenomeNode *gn = (GtGenomeNode*) key;
  GtSpecAspectNodeResult *sanr = (GtSpecAspectNodeResult*) value;
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  const char *tmp;
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

static int gt_spec_results_make_aspect_model(GT_UNUSED void *key, void *value,
                                         void *data, GT_UNUSED GtError *err)
{
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  GtSpecAspect* sa =  (GtSpecAspect*) value;
  GT_UNUSED int rval;
  gt_assert(sa && info);

  info->node_i = 1;
  lua_pushstring(info->L, gt_str_get(sa->name));
  lua_newtable(info->L);
  lua_pushstring(info->L, "successes");
  lua_pushnumber(info->L, sa->successes);
  lua_rawset(info->L, -3);
  lua_pushstring(info->L, "failures");
  lua_pushnumber(info->L, sa->failures);
  lua_rawset(info->L, -3);
  lua_pushstring(info->L, "runtime_errors");
  lua_pushnumber(info->L, sa->runtime_errors);
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

void gt_spec_results_report(GtSpecResults *sr, GtFile *outfile,
                            const char *specfile, bool details, bool colored)
{
  GtSpecResultsReportInfo info;
  gt_assert(sr);

  info.outfile = outfile;
  info.details = details;
  info.colored = colored;

  if (sr->seen_feature || sr->seen_meta || sr->seen_region || sr->seen_comment
        || sr->seen_sequence) {
    gt_file_xprintf(outfile, "According to the specification in %s%s%s,\n\n",
                               colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                               specfile,
                               colored ? GT_SPEC_ANSI_COLOR_RESET : "");
  }
  if (sr->seen_feature) {
    gt_hashmap_foreach_in_key_order(sr->feature_aspects,
                                    gt_spec_results_report_features,
                                    &info, NULL);
  }
  if (sr->seen_meta) {
    gt_file_xprintf(outfile, "a %smeta%s node\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->meta_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  if (sr->seen_region) {
    gt_file_xprintf(outfile, "a %sregion%s node\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->region_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  if (sr->seen_comment) {
    gt_file_xprintf(outfile, "a %scomment%s node\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->comment_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  if (sr->seen_sequence) {
    gt_file_xprintf(outfile, "a %ssequence%s node\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->sequence_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  if (gt_str_array_size(sr->warnings) > 0) {
    gt_file_xprintf(outfile, "\nThere ha%s been %s"GT_WU"%s parser warning%s\n",
                             gt_str_array_size(sr->warnings) > 1 ? "ve" : "s",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             gt_str_array_size(sr->warnings),
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "",
                             gt_str_array_size(sr->warnings) > 1 ? "s." : ".");
    if (details) {
      GtUword i;
      for (i = 0; i < gt_str_array_size(sr->warnings); i++) {
        gt_file_xprintf(outfile, "  - %s\n", gt_str_array_get(sr->warnings, i));
      }

    }
  }
  gt_file_xprintf(outfile, "\nTraversed " GT_WU " CCs "
                           "(" GT_WU " feature types), "
                           "checked " GT_WU " nodes for " GT_WU " aspects.\n",
                           sr->checked_ccs, sr->checked_types,
                           sr->checked_nodes, sr->checked_aspects);
}

/* registry keys for outfile pointer */
static const char *spec_resuserdata = "gt_spec_results_userdata";

static int gt_spec_results_lua_print(lua_State* L) {
    GtFile *outfile;
    int i = 0;

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

int gt_spec_results_render_template(GtSpecResults *sr, const char *template,
                                    GtFile *outfile, const char *specfile,
                                    bool details, bool colored,
                                    const char *runtime_str, GtError *err)
{
  lua_State *L;
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
