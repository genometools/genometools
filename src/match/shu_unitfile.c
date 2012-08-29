/*
  Copyright (c) 2011 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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
#include "lua.h"
#include "lualib.h"

#include "core/assert_api.h"
#include "core/basename_api.h"
#include "core/encseq_api.h"
#include "core/fileutils_api.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/str_api.h"
#include "core/str_array.h"
#include "core/str_array_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"

#include "match/lua_tools.h"

#include "match/shu_unitfile.h"

static int shu_unitfile_load(lua_State *L,
                             const char *unitfile,
                             GtError *err)
{
  int had_err = 0;
  if (luaL_loadfile(L, unitfile) ||
      lua_pcall(L, 0, 0, 0)) {
    had_err = -1;
    gt_error_set(err, "Lua could not load file '%s'!", lua_tostring(L, -1));
  }
  if (!had_err) {
    lua_getglobal(L, "units");
    if (!lua_istable(L, -1)) {
      had_err = -1;
      gt_error_set(err, "Somethings wrong with the unitfile.");
    }
  }
  return had_err;
}

static int shu_unitfile_compare_inner_key_filenames(lua_State *L,
                                                   GtShuUnitFileInfo *unit_info,
                                                   bool *file_set,
                                                   GtStr *mapping_filename,
                                                   unsigned long genome_idx,
                                                   GtError *err)
{
  bool found = false;
  int had_err = 0;
  unsigned long file_idx;
  GtStr *encseq_filename;

  gt_str_reset(mapping_filename);
  gt_str_append_cstr(mapping_filename, lua_tostring(L, -1));

  for (file_idx = 0;
       file_idx < unit_info->num_of_files && !had_err && !found;
       file_idx++) {

    encseq_filename = gt_str_array_get_str(unit_info->file_names, file_idx);
    if (gt_str_cmp(encseq_filename, mapping_filename) == 0) {
      if (file_set[file_idx] != false) {
        gt_error_set(err, "file %s double entry",
                     gt_str_get(mapping_filename));
        had_err = -1;
      }
      else {
        file_set[file_idx] = true;
        unit_info->map_files[file_idx] = genome_idx;
        found = true;
      }
    }
  }
  if (!found && !had_err) {
    gt_error_set(err, "file %s not found in index, part of genome %s",
                 gt_str_get(mapping_filename),
                 gt_str_get(gt_str_array_get_str(unit_info->genome_names,
                                                 genome_idx)));
    had_err = -1;
  }
  return had_err;
}

static int shu_unitfile_traverse_inner_keys(lua_State *L,
                                            GtShuUnitFileInfo *unit_info,
                                            bool *file_set,
                                            GtStr *mapping_filename,
                                            unsigned long genome_idx,
                                            unsigned long *files_added,
                                            GtError *err)
{
  int had_err = 0;

  gt_str_array_add_cstr(unit_info->genome_names, lua_tostring(L, -2));

  lua_pushnil(L); /* the first inner key */
  while (lua_next(L, -2) != 0 && !had_err) {
    had_err = shu_unitfile_compare_inner_key_filenames(L, unit_info, file_set,
                                                       mapping_filename,
                                                       genome_idx,
                                                       err);
    if (!had_err)
      (*files_added)++;
    lua_pop(L, 1);
  }
  return had_err;
}

static int traverse_units(lua_State *L,
                          GtShuUnitFileInfo *unit_info,
                          GtError *err)
{
  int had_err = 0;
  unsigned long genome_idx = 0,
                files_added = 0;
  bool *file_set;
  GtStr *mapping_filename = gt_str_new();

  gt_assert(unit_info->file_names);

  file_set = gt_calloc((size_t) unit_info->num_of_files,
                       sizeof (file_set));

  unit_info->map_files = gt_calloc((size_t) unit_info->num_of_files,
                                   sizeof (unit_info->map_files));

  gt_str_array_reset(unit_info->genome_names);

  lua_pushnil(L); /*the first outer key*/
  while (lua_next(L, -2) != 0 && !had_err) {
    had_err = shu_unitfile_traverse_inner_keys(L, unit_info, file_set,
                                               mapping_filename,
                                               genome_idx, &files_added, err);
    lua_pop(L, 1);
    genome_idx++;
  }
  lua_pop(L, 1);

  gt_free(file_set);

  if (!had_err) {
    unit_info->num_of_genomes = genome_idx;
    gt_assert(unit_info->num_of_genomes ==
        gt_str_array_size(unit_info->genome_names));
  }
  if (!had_err && files_added != unit_info->num_of_files) {
    had_err = -1;
    gt_error_set(err, "number of files in index (%lu) and unitfile (%lu)! "
                 "differ!", unit_info->num_of_files, files_added);
  }
  if (!had_err) {
    unsigned long file_idx;
    for (file_idx = 0; file_idx < unit_info->num_of_files; file_idx++) {
      gt_log_log("file: %lu belongs to genome: %s",
                 file_idx,
                 gt_str_array_get(unit_info->genome_names,
                                  unit_info->map_files[file_idx]));
    }
  }
  gt_str_delete(mapping_filename);
  return had_err;
}

int gt_shu_unit_file_info_read(const GtStr *unitfile,
                               GtShuUnitFileInfo *unit_info,
                               GT_UNUSED GtLogger *logger,
                               GtError *err)
{
  int had_err = 0;
  /* open Lua */
  lua_State *L = luaL_newstate();
  /* library needed?
  lua_openlibs(L); */

  had_err = shu_unitfile_load(L, gt_str_get(unitfile), err);
  if (!had_err)
    had_err = traverse_units(L, unit_info, err);

  lua_close(L);
  return had_err;
}

void gt_shu_unit_info_delete(GtShuUnitFileInfo *unit_info) {
  if (unit_info != NULL) {
    gt_free(unit_info->map_files);
    gt_str_array_delete(unit_info->genome_names);
    gt_free(unit_info);
  }
}

static void gt_shu_unit_info_files_as_units(GtShuUnitFileInfo *unit_info)
{
  unsigned long i_idx;

  unit_info->num_of_genomes = unit_info->num_of_files;
  unit_info->genome_names = gt_str_array_new();
  for (i_idx = 0; i_idx < unit_info->num_of_files; i_idx++)
  {
    gt_str_array_add(unit_info->genome_names,
                     gt_str_array_get_str(unit_info->file_names, i_idx));
  }
}

GtShuUnitFileInfo *gt_shu_unit_info_new(const GtEncseq *encseq)
{
  GtShuUnitFileInfo *unit_info = gt_malloc(sizeof (*unit_info));

  unit_info->map_files = NULL;
  unit_info->genome_names = NULL;
  unit_info->num_of_genomes = 0;
  unit_info->num_of_files = gt_encseq_num_of_files(encseq);
  unit_info->file_names = gt_encseq_filenames(encseq);
  unit_info->encseq = encseq;
  gt_shu_unit_info_files_as_units(unit_info);
  return unit_info;
}
