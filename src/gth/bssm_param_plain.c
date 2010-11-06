/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "core/ma_api.h"
#include "core/str_api.h"
#include "core/xansi_api.h"
#include "gth/bssm_param_plain.h"
#include "gth/bssm_param_rep.h"

#define PRECISION  7

static void write_model(GtStr *str, const char *model_cstr,
                        const GthBSSMModel *model)
{
  unsigned long i, j, k, l;
  gt_assert(str && model_cstr && model);
  gt_assert(model->hypothesis_num == 2 || model->hypothesis_num == 7);
  gt_str_append_cstr(str, "  ");
  gt_str_append_cstr(str, model_cstr);
  gt_str_append_cstr(str, " = {\n");
  gt_str_append_cstr(str, "    hypothesis_num = ");
  gt_str_append_ulong(str, model->hypothesis_num);
  gt_str_append_cstr(str, ",\n");
  gt_str_append_cstr(str, "    window_size_left = ");
  gt_str_append_ulong(str, model->window_size_left);
  gt_str_append_cstr(str, ",\n");
  gt_str_append_cstr(str, "    window_size_right = ");
  gt_str_append_ulong(str, model->window_size_right);
  gt_str_append_cstr(str, ",\n");
  for (i = 0; i < model->hypothesis_num; i++) {
    gt_str_append_cstr(str, "    {\n");
    for (j = 0; j < WINSIZE + 2; j++) {
      gt_str_append_cstr(str, "      {\n");
      for (k = 0; k < 4; k++) {
        gt_str_append_cstr(str, "        { ");
        for (l = 0; l < 4; l++) {
          if (l)
            gt_str_append_cstr(str, ", ");
          if (model->hypothesis_num == 2) {
            gt_str_append_double(str, model->hypotables.hypo2table[i][j][k][l],
                                 PRECISION);
          }
          else {
            gt_str_append_double(str, model->hypotables.hypo7table[i][j][k][l],
                                 PRECISION);
          }
        }
        gt_str_append_cstr(str, " },\n");
      }
      gt_str_append_cstr(str, "      },\n");
    }
    gt_str_append_cstr(str, "    },\n");
  }
  gt_str_append_cstr(str, "  }");
}

void gth_bssm_param_plain_write(const GthBSSMParam *bssm_param, FILE *outfp)
{
  GtStr *str;
  gt_assert(bssm_param && outfp);
  str = gt_str_new();
  gt_str_append_cstr(str, "BSSM = {\n");
  if (bssm_param->gt_donor_model_set) {
    write_model(str, "gt_donor_model", &bssm_param->gt_donor_model);
    gt_str_append_cstr(str, ",\n");
  }
  if (bssm_param->gc_donor_model_set) {
    write_model(str, "gc_donor_model", &bssm_param->gc_donor_model);
    gt_str_append_cstr(str, ",\n");
  }
  if (bssm_param->ag_acceptor_model_set) {
    write_model(str, "ag_acceptor_model", &bssm_param->ag_acceptor_model);
    gt_str_append_char(str, '\n');
  }
  gt_str_append_cstr(str, "}\n");
  gt_xfwrite(gt_str_get(str), sizeof (char), gt_str_length(str), outfp);
  gt_str_delete(str);
}

static int read_bssm_model(GthBSSMModel *bssm_model, lua_State *L, GtError *err)
{
#ifndef NDEBUG
  int stack_size;
#endif
  unsigned long i, j, k, l;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(bssm_model && L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  /* read hypothesis number */
  lua_pushliteral(L, "hypothesis_num");
  lua_gettable(L, -2);
  if (!lua_isnumber(L, -1)) {
    gt_error_set(err, "'hypothesis_num' could not be found in BSSM file or is "
                      "not a number");
    had_err = -1;
    lua_pop(L, 1);
  }
  if (!had_err) {
    bssm_model->hypothesis_num = lua_tonumber(L, -1);
    lua_pop(L, 1);
    if (!(bssm_model->hypothesis_num == 2 || bssm_model->hypothesis_num == 7)) {
      gt_error_set(err, "'hypothesis_num' in BSSM file must be either 2 or 7");
      had_err = -1;
    }
  }
  /* read window size left */
  if (!had_err) {
    lua_pushliteral(L, "window_size_left");
    lua_gettable(L, -2);
    if (!lua_isnumber(L, -1)) {
      gt_error_set(err, "'window_size_left' could not be found in BSSM file or "
                        "is not a number");
      had_err = -1;
      lua_pop(L, 1);
    }
    if (!had_err) {
      bssm_model->window_size_left = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }
  }
  /* read window size right */
  if (!had_err) {
    lua_pushliteral(L, "window_size_right");
    lua_gettable(L, -2);
    if (!lua_isnumber(L, -1)) {
      gt_error_set(err, "'window_size_right' could not be found in BSSM file "
                        "or is not a number");
      had_err = -1;
      lua_pop(L, 1);
    }
    if (!had_err) {
      bssm_model->window_size_right = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }
  }
  /* read hypothesis table */
  for (i = 0; !had_err && i < bssm_model->hypothesis_num; i++) {
    lua_pushnumber(L, i+1);
    lua_gettable(L, -2);
    if (!lua_istable(L, -1)) {
      gt_error_set(err, "incomplete BSSM");
      had_err = -1;
    }
    for (j = 0; !had_err && j < WINSIZE + 2; j++) {
      lua_pushnumber(L, j+1);
      lua_gettable(L, -2);
      if (!lua_istable(L, -1)) {
        gt_error_set(err, "incomplete BSSM");
        had_err = -1;
      }
      for (k = 0; !had_err && k < 4; k++) {
        lua_pushnumber(L, k+1);
        lua_gettable(L, -2);
        if (!lua_istable(L, -1)) {
          gt_error_set(err, "incomplete BSSM");
          had_err = -1;
        }
        for (l = 0; !had_err && l < 4; l++) {
          lua_pushnumber(L, l+1);
          lua_gettable(L, -2);
          if (!lua_isnumber(L, -1)) {
            gt_error_set(err, "incomplete BSSM");
            had_err = -1;
          }
          else if (bssm_model->hypothesis_num == 2)
            bssm_model->hypotables.hypo2table[i][j][k][l] = lua_tonumber(L, -1);
          else
            bssm_model->hypotables.hypo7table[i][j][k][l] = lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
        lua_pop(L, 1);
      }
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
  }
  gt_assert(lua_gettop(L) == stack_size);
  return had_err;
}

static int fill_bssm_param(GthBSSMParam *bssm_param, lua_State *L, GtError *err)
{
#ifndef NDEBUG
  int stack_size;
#endif
  bool model_read = false;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(bssm_param && L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  lua_getglobal(L, "BSSM");
  gt_assert(lua_istable(L, -1));
  bssm_param->version_num = BSSMPARAMVERSION;
  /* try to read gt donor model */
  lua_pushliteral(L, "gt_donor_model");
  lua_gettable(L, -2);
  if (lua_istable(L, -1)) {
    had_err = read_bssm_model(&bssm_param->gt_donor_model, L, err);
    if (!had_err) {
      bssm_param->gt_donor_model_set = true;
      model_read = true;
    }
  }
  lua_pop(L, 1);
  /* try to read gc donor model */
  if (!had_err) {
    lua_pushliteral(L, "gc_donor_model");
    lua_gettable(L, -2);
    if (lua_istable(L, -1)) {
      had_err = read_bssm_model(&bssm_param->gc_donor_model, L, err);
      if (!had_err) {
        bssm_param->gc_donor_model_set = true;
        model_read = true;
      }
    }
    lua_pop(L, 1);
  }
  /* try to read ag acceptor model */
  if (!had_err) {
    lua_pushliteral(L, "ag_acceptor_model");
    lua_gettable(L, -2);
    if (lua_istable(L, -1)) {
      had_err = read_bssm_model(&bssm_param->ag_acceptor_model, L, err);
      if (!had_err) {
        bssm_param->ag_acceptor_model_set = true;
        model_read = true;
      }
    }
    lua_pop(L, 1);
  }
  /* make sure at least one model has been read */
  if (!had_err && !model_read) {
    gt_error_set(err, "BSSM file contains no model");
    had_err = -1;
  }
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return had_err;
}

static int read_bssm_file(GthBSSMParam *bssm_param, lua_State *L,
                          const char *filename, GtError *err)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int had_err = 0;
  gt_error_check(err);
  gt_assert(L && filename);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  if (luaL_loadfile(L, filename) || lua_pcall(L, 0, 0, 0)) {
    gt_error_set(err, "cannot read plain BSSM file: %s", lua_tostring(L, -1));
    had_err = -1;
    lua_pop(L, 1);
  }
  if (!had_err) {
    lua_getglobal(L, "BSSM");
    if (lua_isnil(L, -1) || !lua_istable(L, -1)) {
      gt_error_set(err, "'BSSM' is not defined or is not a table in \"%s\"",
                   filename);
      had_err = -1;
    }
    lua_pop(L, 1);
  }
  if (!had_err)
    had_err = fill_bssm_param(bssm_param, L, err);
  gt_assert(lua_gettop(L) == stack_size);
  return had_err;
}

GthBSSMParam* gth_bssm_param_plain_read(const char *filename, GtError *err)
{
  GthBSSMParam *bssm_param;
  lua_State *L;
  int had_err;
  gt_error_check(err);
  gt_assert(filename);
  bssm_param = gt_calloc(1, sizeof *bssm_param);
  L = luaL_newstate();
  had_err = read_bssm_file(bssm_param, L, filename, err);
  lua_close(L);
  if (had_err) {
    gt_free(bssm_param);
    return NULL;
  }
  return bssm_param;
}
