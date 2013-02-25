/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2005 Michael E Sparks <mespar1@iastate.edu>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/bioseq.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/ma_api.h"
#include "core/xansi_api.h"
#include "gth/bssm_param_hard_coded.h"
#include "gth/bssm_param.h"
#include "gth/gthoutput.h"
#include "gth/gthprobdef.h"
#include "gth/gthspeciestab.h"
#include "gth/showbool.h"

/*
  This is a collection of functions associated with
  manipulating bssm_param objects for training purposes.

  Input data files are--strictly!--named as follows, using
  an obvious schema:
    F0_don  F1_don  F2_don  Fi_don  T0_don  T1_don  T2_don
    F0_acc  F1_acc  F2_acc  Fi_acc  T0_acc  T1_acc  T2_acc
  Phase is denoted as follows:
    1 -> C O D |
    2 -> C | O D
    0 -> C O | D
*/

#define PSEUDOPROB      0.05
#define NULLPROB        0.0
#define INITVAL_INT     0
#define MYVERSION       2  /* Version of BSSM param. Check with
                              BSSMPARAMVERSION in bssm_param.h */

#define NUMOFFILES (sizeof (filenames)/sizeof (filenames[0]))

#define BSSMENVNAME  "BSSMDIR"

#define BSSM_PRECISION  8

static char *filenames[] =
{
  "T1",
  "T2",
  "T0",
  "F1",
  "F2",
  "F0",
  "Fi"
};

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
                                 BSSM_PRECISION);
          }
          else {
            gt_str_append_double(str, model->hypotables.hypo7table[i][j][k][l],
                                 BSSM_PRECISION);
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

static void bssm_param_plain_write(const GthBSSMParam *bssm_param, FILE *outfp)
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

static GthBSSMParam* bssm_param_plain_read(const char *filename, GtError *err)
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

static int bssm_model_read(GthBSSMModel *bssm_model, FILE *file, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_xfread(&bssm_model->hypothesis_num, sizeof (unsigned long), 1, file);
  if (bssm_model->hypothesis_num != HYPOTHESIS7 &&
      bssm_model->hypothesis_num != HYPOTHESIS2) {
    gt_error_set(err, "BSSM model contains unknown hypothesis number");
    had_err = -1;
  }
  if (!had_err) {
    gt_xfread(&bssm_model->window_size_left, sizeof (unsigned long), 1, file);
    gt_xfread(&bssm_model->window_size_right, sizeof (unsigned long), 1, file);
    switch (bssm_model->hypothesis_num) {
      case HYPOTHESIS2:
        gt_xfread(&bssm_model->hypotables.hypo2table, sizeof (Hypo2table), 1,
                  file);
        break;
      case HYPOTHESIS7:
        gt_xfread(&bssm_model->hypotables.hypo7table, sizeof (Hypo7table), 1,
                  file);
        break;
      default: gt_assert(0);
    }
  }
  return had_err;
}

GthBSSMParam* gth_bssm_param_new(void)
{
  return gt_calloc(1, sizeof (GthBSSMParam));
}

static GthBSSMParam* load_old_binary_format(GtStr *path, const char *filename,
                                            GtError *err)
{
  GthBSSMParam *bssm_param;
  int had_err = 0;
  FILE *file;
  gt_error_check(err);
  gt_assert(path && filename);

  file = gt_xfopen(gt_str_get(path), "r");

  /* read version number and check if equals version number 2 */
  bssm_param = gt_malloc(sizeof *bssm_param);
  gt_xfread(&bssm_param->version_num,  sizeof (unsigned char), 1, file);
  if (bssm_param->version_num != (unsigned char) 2) {
    gt_error_set(err, "BSSM file %s has unrecognized version number %u",
                 filename, bssm_param->version_num);
    had_err = -1;
  }

  if (!had_err) {
    /* read in model variables */
    gt_xfread(&bssm_param->gt_donor_model_set,  sizeof (bool), 1, file);
    gt_xfread(&bssm_param->gc_donor_model_set,  sizeof (bool), 1, file);
    gt_xfread(&bssm_param->ag_acceptor_model_set,  sizeof (bool), 1, file);

    /* check if at least one model is set in file */
    if (!bssm_param->gt_donor_model_set &&
        !bssm_param->gc_donor_model_set &&
        !bssm_param->ag_acceptor_model_set) {
      gt_error_set(err, "BSSM file %s apparently contains no model", filename);
      had_err = -1;
    }
  }

  /* read GT donor site model */
  if (!had_err && bssm_param->gt_donor_model_set)
    had_err = bssm_model_read(&bssm_param->gt_donor_model, file, err);

  /* read GC donor site model */
  if (!had_err && bssm_param->gc_donor_model_set)
    had_err = bssm_model_read(&bssm_param->gc_donor_model, file, err);

  /* read AG acceptor site model */
  if (!had_err && bssm_param->ag_acceptor_model_set)
    had_err = bssm_model_read(&bssm_param->ag_acceptor_model, file, err);

  gt_xfclose(file);

  if (had_err) {
    gth_bssm_param_delete(bssm_param);
    return NULL;
  }
  return bssm_param;
}

GthBSSMParam* gth_bssm_param_load(const char *filename, GtError *err)
{
  GthBSSMParam *bssm_param = NULL;
  GtStr *path = gt_str_new();
  int had_err = 0;

  gt_error_check(err);

  if (gt_file_exists(filename))
    gt_str_append_cstr(path, filename);
  else {
    if (strchr(filename, '/')) {
      gt_error_set(err, "filename \"%s\" contains illegal symbol '/': the path "
                        "list specified by environment variable \"%s\" cannot "
                        "be searched for it", filename, BSSMENVNAME);
      had_err = -1;
    }
    if (!had_err)
      had_err = gt_file_find_in_env(path, filename, BSSMENVNAME, err);
    if (!had_err && !gt_str_length(path)) {
      gt_error_set(err, "file \"%s\" not found in directory list specified by "
                        "environment variable %s", filename, BSSMENVNAME);
      had_err = -1;
    }
    if (!had_err) {
      /* path found -> append filename */
      gt_str_append_char(path, '/');
      gt_str_append_cstr(path, filename);
    }
  }

  if (!had_err) {
    if (!(bssm_param = bssm_param_plain_read(gt_str_get(path), err)))
      had_err = -1;
    if (had_err) {
      /* loading new plain text format didn't work -> try old binary format */
      if ((bssm_param = load_old_binary_format(path, filename, NULL))) {
        /* loading binary format worked -> unset error */
        gt_error_unset(err);
        had_err = 0;
      }
    }
  }

  gt_str_delete(path);

  if (had_err) {
    gth_bssm_param_delete(bssm_param);
    return NULL;
  }
  return bssm_param;
}

GthBSSMParam* gth_bssm_param_extract(unsigned long speciesnum, GtError *err)
{
  GthBSSMParam *bssm_param;
  unsigned long i, j, k, l;

  gt_error_check(err);

  bssm_param = gth_bssm_param_new();
  bssm_param->version_num = (unsigned char) 2;
  bssm_param->gt_donor_model_set    = true;
  bssm_param->gc_donor_model_set    = false;
  bssm_param->ag_acceptor_model_set = true;

  if (speciesnum <= 1) {
    /* read in the human and mouse cases */
    bssm_param->gt_donor_model.hypothesis_num    = HYPOTHESIS2;
    bssm_param->ag_acceptor_model.hypothesis_num = HYPOTHESIS2;
    for (i = 0; i < HYPOTHESIS2; i++) {
      for (j = 0; j < WINSIZE + 2; j++) {
        for (k = 0; k < 4; k++) {
          for (l = 0; l < 4; l++) {
            bssm_param->gt_donor_model.hypotables.hypo2table[i][j][k][l] =
              (GthFlt) GU_2[speciesnum][i][j][k][l];
            bssm_param->ag_acceptor_model.hypotables.hypo2table[i][j][k][l] =
              (GthFlt) AG_2[speciesnum][i][j][k][l];
          }
        }
      }
    }
  }
  else if (speciesnum <  NUMOFSPECIES) {
    /* read in all others */
    bssm_param->gt_donor_model.hypothesis_num    = HYPOTHESIS7;
    bssm_param->ag_acceptor_model.hypothesis_num = HYPOTHESIS7;
    for (i = 0; i < HYPOTHESIS7; i++) {
      for (j = 0; j < WINSIZE + 2; j++) {
        for (k = 0; k < 4; k++) {
          for (l = 0; l < 4; l++) {
            bssm_param->gt_donor_model.hypotables.hypo7table[i][j][k][l] =
              (GthFlt) GU_7[speciesnum - 2][i][j][k][l];
            bssm_param->ag_acceptor_model.hypotables.hypo7table[i][j][k][l] =
              (GthFlt) AG_7[speciesnum - 2][i][j][k][l];
          }
        }
      }
    }
  }
  else {
    gt_error_set(err, "illegal speciesnum (speciesnum == %lu)", speciesnum);
    gth_bssm_param_delete(bssm_param);
    return NULL;
  }

  /* setting the window sizes */
  bssm_param->gt_donor_model.window_size_left     = wsize[speciesnum][0][0];
  bssm_param->gt_donor_model.window_size_right    = wsize[speciesnum][0][1];
  bssm_param->ag_acceptor_model.window_size_left  = wsize[speciesnum][1][0];
  bssm_param->ag_acceptor_model.window_size_right = wsize[speciesnum][1][1];

  return bssm_param;
}

void gth_bssm_param_delete(GthBSSMParam *bssm_param)
{
  if (!bssm_param) return;
  gt_free(bssm_param);
}

int gth_bssm_param_save(GthBSSMParam *bssm_param, const char *filename,
                        GtError *err)
{
  FILE *file;
  int had_err = 0;

  gt_error_check(err);

  file = gt_fa_xfopen(filename, "w");
  gt_assert(file);

  /* check if at least one model is set */
  if (!bssm_param->gt_donor_model_set &&
      !bssm_param->gc_donor_model_set &&
      !bssm_param->ag_acceptor_model_set) {
    gt_error_set(err, "BSSM parameter to write contain no model");
    had_err = -1;
  }

  if (!had_err)
    bssm_param_plain_write(bssm_param, file);

  gt_fa_xfclose(file);

  return had_err;
}

static bool bssm_model_is_seven_class(const GthBSSMModel *bssm_model)
{
  gt_assert(bssm_model);
  return bssm_model->hypothesis_num == HYPOTHESIS7;
}

bool gth_bssm_param_is_seven_class(const GthBSSMParam  *bssm_param)
{
  gt_assert(bssm_param);
  return (!bssm_param->gt_donor_model_set ||
          bssm_model_is_seven_class(&bssm_param->gt_donor_model)) &&
         (!bssm_param->gc_donor_model_set ||
          bssm_model_is_seven_class(&bssm_param->gc_donor_model)) &&
         (!bssm_param->ag_acceptor_model_set ||
          bssm_model_is_seven_class(&bssm_param->ag_acceptor_model));
}

/* The following function outouts <bssm_model>.
   It is assumed that the Hypo7table is used. */
static void bssm_model_echo(const GthBSSMModel *bssm_model, FILE *outfp)
{
  unsigned long i, j, k, l;

  gt_assert(bssm_model_is_seven_class(bssm_model));

  for (i = 0; i < HYPOTHESIS7; i++) {
    fprintf(outfp,"\n\nHypothesis: %lu", i);
    for (j = 0; j < STRINGSIZE; j++) {
      fprintf(outfp,"\n");
      for (k = 0; k < ALPHSIZE; k++) {
        fprintf(outfp,"\n");
        for (l = 0; l < ALPHSIZE; l++) {
          fprintf(outfp,"%.4f ", bssm_model->hypotables.hypo7table[i][j][k][l]);
        }
      }
    }
  }
  fprintf(outfp,"\n\n");
}

void gth_bssm_param_echo(const GthBSSMParam *bssm_param, FILE *outfp)
{
  gt_assert(bssm_param && outfp);
  fprintf(outfp,"BSSMPARAMVERSION is %u\n\n", bssm_param->version_num);
  fprintf(outfp,"Is the GT donor model set? -> %s\n",
          GTH_SHOWBOOL(bssm_param->gt_donor_model_set));
  fprintf(outfp,"Is the GC donor model set? -> %s\n\n",
          GTH_SHOWBOOL(bssm_param->gc_donor_model_set));
  fprintf(outfp,"Is the AG acceptor model set? -> %s\n\n",
          GTH_SHOWBOOL(bssm_param->ag_acceptor_model_set));

  if (bssm_param->gt_donor_model_set) {
    fprintf(outfp,"reporting GT donor model parameterization");
    bssm_model_echo(&bssm_param->gt_donor_model, outfp);
  }

  if (bssm_param->gc_donor_model_set) {
    fprintf(outfp,"reporting GC donor model parameterization");
    bssm_model_echo(&bssm_param->gc_donor_model, outfp);
  }

  if (bssm_param->ag_acceptor_model_set) {
    fprintf(outfp,"reporting AG acceptor model parameterization");
    bssm_model_echo(&bssm_param->ag_acceptor_model, outfp);
  }
}

void gth_bssm_param_show_info(const GthBSSMParam *bssm_param, GtFile *outfp)
{
#define SEVENCLASSSTRING        "seven-class"
#define TWOCLASSSTRING          "two-class"

#define PRINT_CLASS_STRING(MODEL) \
  if (bssm_param->MODEL##_model_set) \
  { \
    gt_file_xprintf(outfp, " (%s)", \
                    bssm_param->MODEL##_model.hypothesis_num == HYPOTHESIS7 \
                    ? SEVENCLASSSTRING : TWOCLASSSTRING); \
  } \
  gt_file_xfputc('\n', outfp);

  gt_file_xprintf(outfp,
                  "%c the specified BSSM parameter file contains the following "
                  "models:\n", COMMENTCHAR);
  gt_file_xprintf(outfp, "%c GT donor sites   = %s", COMMENTCHAR,
                  GTH_SHOWBOOL(bssm_param->gt_donor_model_set));
  PRINT_CLASS_STRING(gt_donor);

  gt_file_xprintf(outfp, "%c GC donor sites   = %s", COMMENTCHAR,
                  GTH_SHOWBOOL(bssm_param->gc_donor_model_set));
  PRINT_CLASS_STRING(gc_donor);

  gt_file_xprintf(outfp, "%c AG acceptor sites= %s", COMMENTCHAR,
                  GTH_SHOWBOOL(bssm_param->ag_acceptor_model_set));
  PRINT_CLASS_STRING(ag_acceptor);
}

static void set_window_sizes_in_Bssmmodel(GthBSSMModel *bssm_model)
{
  bssm_model->hypothesis_num = HYPOTHESIS7;

  /* We have decided to leave maximum splice signal
     extent window at...a maximal value!            */
  bssm_model->window_size_left  = MAXSPLICESIG;
  bssm_model->window_size_right = MAXSPLICESIG;
}

/* updates the BSSM parameterization file */
static void build_bssm(GtBioseq *bioseq, GthBSSMModel *bssm_model,
                       unsigned int hypothesisnum)
{
  unsigned long mono_ct[STRINGSIZE-1][ALPHSIZE],         /* Mononuc freq */
                di_ct[STRINGSIZE-1][ALPHSIZE][ALPHSIZE]; /* Dinuc freq */
  double mono_freq,      /* Mononuc relative freq */
         di_freq;        /* Dinuc relative freq */
  unsigned long i, j, k, /* Iterator variables */
                len, curlen = 0,
                num_entries = gt_bioseq_number_of_sequences(bioseq);
  GtUchar *encoded_seq = NULL;

  /* Inits of local variables */
  for (i = 0; i < (STRINGSIZE-1); i++) {
    for (j = 0; j < ALPHSIZE; j++) {
      mono_ct[i][j] = INITVAL_INT;
      for (k = 0; k < ALPHSIZE; k++)
        di_ct[i][j][k] = INITVAL_INT;
    }
  }

  /* mononucleotides */
  for (j = 0; j < num_entries; j++) {
    len = gt_bioseq_get_sequence_length(bioseq, j);
    gt_assert(len == STRINGSIZE);
    if (len > curlen) {
      encoded_seq = gt_realloc(encoded_seq, len);
      curlen = len;
    }
    gt_bioseq_get_encoded_sequence(bioseq, encoded_seq, j);
    for (i = 0; i < (STRINGSIZE-1); i++) {
      gt_assert(encoded_seq[i] < ALPHSIZE);
      mono_ct[i][encoded_seq[i]]++;
    }
  }

  /* dinucleotides */
  for (j = 0; j < num_entries; j++) {
    len = gt_bioseq_get_sequence_length(bioseq, j);
    gt_assert(len == STRINGSIZE);
    if (len > curlen) {
      encoded_seq = gt_realloc(encoded_seq, len);
      curlen = len;
    }
    gt_bioseq_get_encoded_sequence(bioseq, encoded_seq, j);
    for (i = 0; i < (STRINGSIZE-1); i++) {
      di_ct[i][encoded_seq[i]]
              [encoded_seq[i + 1]]++;
    }
  }

  gt_free(encoded_seq);

  /* Record equilibrium frequencies (1st ``slot" in transition freqs) */
  for (i = 0; i < ALPHSIZE; i++) {
    for (j = 0; j < ALPHSIZE; j++) {
      bssm_model->hypotables
      .hypo7table[hypothesisnum][0][i][j] = (GthFlt)
                                            mono_ct[0][i] / num_entries;
    }
  }

  /* Populate the remaining transition frequencies */
  for (k = 1; k < STRINGSIZE; k++) {
    for (i = 0; i < ALPHSIZE; i++) {
      mono_freq = (double) mono_ct[k-1][i] / num_entries;
      for (j = 0; j < ALPHSIZE; j++) {
        di_freq = (double) di_ct[k-1][i][j] / num_entries;
        if (mono_freq == 0.0) {
          bssm_model->hypotables
          .hypo7table[hypothesisnum][k][i][j] = (GthFlt) NULLPROB;
        }
        else {
          bssm_model->hypotables
          .hypo7table[hypothesisnum][k][i][j] = (GthFlt)
                                                (di_freq / mono_freq);
        }
      }

      /* Remove non-zero transition probabilities:
         Briefly, 0.0 entries (dinucleotide absent in training corpus) are
         replaced arbitrarily by PSEUDOPROB, and non-0.0 entries p are replaced
         by p = p * (1 - 4 * PSEUDOPROB) + PSEUDOPROB */
      for (j = 0; j < ALPHSIZE; ++j) {
        /* If any entry is NULLPROB, ALL elements in the row need fixed */
        if (bssm_model->hypotables
            .hypo7table[hypothesisnum][k][i][j] == NULLPROB) {
          /* Fix all elements in the row, then break */
          for (j = 0; j < ALPHSIZE; j++) {
            if (bssm_model->hypotables
                .hypo7table[hypothesisnum][k][i][j] == NULLPROB) {
               bssm_model->hypotables
                .hypo7table[hypothesisnum][k][i][j] = (GthFlt)
                                                      PSEUDOPROB;
            }
            else {
              /* Adjust non-zero transition prob */
              bssm_model->hypotables.hypo7table[hypothesisnum][k][i][j] =
                (GthFlt)
                (bssm_model->hypotables.hypo7table[hypothesisnum][k][i][j] *
                 (1 - (4 * PSEUDOPROB)) + PSEUDOPROB);
            }
          }
          break;
        }
      }
    }
  }
}

int gth_bssm_param_parameterize(GthBSSMParam *bssm_param, const char *path,
                                Termtype termtype, bool gzip, GtError *err)
{
  GtAlphabet *alphabet = NULL;
  GtBioseq *bioseq;
  GtStr *file2proc;
  unsigned long i, j;
  int had_err = 0;
  gt_error_check(err);

  file2proc = gt_str_new();

  /* set version number */
  bssm_param->version_num = (unsigned char) MYVERSION;

  /* set model to true and set window sizes */
  switch (termtype) {
    case GT_DONOR_TYPE:
      bssm_param->gt_donor_model_set = true;
      set_window_sizes_in_Bssmmodel(&bssm_param->gt_donor_model);
      break;
    case GC_DONOR_TYPE:
      bssm_param->gc_donor_model_set = true;
      set_window_sizes_in_Bssmmodel(&bssm_param->gc_donor_model);
      break;
    case AG_ACCEPTOR_TYPE:
      bssm_param->ag_acceptor_model_set = true;
      set_window_sizes_in_Bssmmodel(&bssm_param->ag_acceptor_model);
      break;
    default: gt_assert(0);
  }

  for (i = 0; !had_err && i < NUMOFFILES; i++) {
    /* process datafile */
    gt_str_append_cstr(file2proc, path);
    switch (termtype) {
      case GT_DONOR_TYPE:
        gt_str_append_cstr(file2proc, "/GT_donor/");
        gt_str_append_cstr(file2proc, filenames[i]);
        break;
      case GC_DONOR_TYPE:
        gt_str_append_cstr(file2proc, "/GC_donor/");
        gt_str_append_cstr(file2proc, filenames[i]);
        break;
      case AG_ACCEPTOR_TYPE:
        gt_str_append_cstr(file2proc, "/AG_acceptor/");
        gt_str_append_cstr(file2proc, filenames[i]);
        break;
      default: gt_assert(0);
    }

    if (gzip)
      gt_str_append_cstr(file2proc, ".gz");

    if (!(bioseq = gt_bioseq_new(gt_str_get(file2proc), err)))
      had_err = -1;

    if (!had_err)
      alphabet = gt_bioseq_get_alphabet(bioseq);

    /* check here if all sequences have the length 102 and correct bases at
       positions 51 and 52 (i.e., GT, GC, or AG) */
    for (j = 0; !had_err && j < gt_bioseq_number_of_sequences(bioseq); j++) {
      GtUchar encoded_seq[2];
      /* check length */
      if (gt_bioseq_get_sequence_length(bioseq, j) != STRINGSIZE) {
        gt_error_set(err, "sequence %lu in file \"%s\" does not have length %u",
                     j, gt_str_get(file2proc), STRINGSIZE);
        had_err = -1;
      }
      encoded_seq[0] = gt_bioseq_get_encoded_char(bioseq, j, 50);
      encoded_seq[1] = gt_bioseq_get_encoded_char(bioseq, j, 51);
      if (!had_err) {
        /* check base correctness */
        switch (termtype) {
          case GT_DONOR_TYPE:
            if (encoded_seq[0] != gt_alphabet_encode(alphabet, 'G') ||
                encoded_seq[1] != gt_alphabet_encode(alphabet, 'T')) {
              gt_error_set(err, "sequence %lu in file \"%s\" is not a GT "
                                "sequence", j, gt_str_get(file2proc));
              had_err = -1;
            }
            break;
          case GC_DONOR_TYPE:
            if (encoded_seq[0] != gt_alphabet_encode(alphabet, 'G') ||
                encoded_seq[1] != gt_alphabet_encode(alphabet, 'C')) {
              gt_error_set(err, "sequence %lu in file \"%s\" is not a GC "
                                "sequence", j, gt_str_get(file2proc));
              had_err = -1;
            }
            break;
          case AG_ACCEPTOR_TYPE:
            if (encoded_seq[0] != gt_alphabet_encode(alphabet, 'A') ||
                encoded_seq[1] != gt_alphabet_encode(alphabet, 'G')) {
              gt_error_set(err, "sequence %lu in file \"%s\" is not a AG "
                                "sequence", j, gt_str_get(file2proc));
              had_err = -1;
            }
            break;
          default: gt_assert(0);
        }
      }
    }

    if (!had_err) {
      switch (termtype) {
        case GT_DONOR_TYPE:
          build_bssm(bioseq, &bssm_param->gt_donor_model, i);
          break;
        case GC_DONOR_TYPE:
          build_bssm(bioseq, &bssm_param->gc_donor_model, i);
          break;
        case AG_ACCEPTOR_TYPE:
          build_bssm(bioseq, &bssm_param->ag_acceptor_model, i);
          break;
        default: gt_assert(0);
      }
    }

    /* reset */
    gt_str_reset(file2proc);

    /* free space */
    gt_bioseq_delete(bioseq);
  }
  gt_str_delete(file2proc);

  return had_err;
}
