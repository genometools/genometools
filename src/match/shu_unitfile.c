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
#include "core/ma_api.h"
#include "core/str_api.h"
#include "core/str_array.h"
#include "core/str_array_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"

#include "match/lua_tools.h"

#include "match/shu_unitfile.h"

static int load(lua_State *L,
                 const char *unitfile,
                 GtError *err)
{
  int had_err = 0;
  if (luaL_loadfile(L, unitfile) ||
      lua_pcall(L, 0, 0, 0))
  {
    had_err = 1;
    gt_error_set(err, "Lua could not load file '%s'!", lua_tostring(L, -1));
  }
  if (!had_err)
  {
    lua_getglobal(L, "units");
    if (!lua_istable(L, -1))
    {
      had_err = 1;
      gt_error_set(err, "Somethings wrong with the unitfile.");
    }
  }
  return had_err;
}

static int traverse_units(lua_State *L,
                          struct GtShuUnitFileInfo_tag *unit_info,
                          GtError *err)
{
  int had_err = 0,
      found = 0;
  unsigned long file_idx,
                genome_idx = 0,
                files_added = 0;
  bool *file_set;

  gt_assert(unit_info->file_names);

  file_set = gt_calloc((size_t) unit_info->num_of_files,
                       sizeof (file_set));

  unit_info->map_files = gt_calloc((size_t) unit_info->num_of_files,
                                   sizeof (unit_info->map_files));
  unit_info->genome_names = gt_str_array_new();

  lua_pushnil(L); /*the first outer key*/
  while (lua_next(L, -2) != 0 && !had_err)
  {
    gt_str_array_add_cstr(unit_info->genome_names, lua_tostring(L, -2));
    /*fprintf(stderr, "%s\n", gt_str_array_get(unit_info->genome_names,
                                             genome_idx));*/

    lua_pushnil(L); /* the first inner key */
    while (lua_next(L, -2) != 0 && !had_err)
    {
      GtStr *mapping_filename;

      mapping_filename = gt_str_new_cstr(lua_tostring(L, -1));
      /*fprintf(stderr, "%s\n", gt_str_get(mapping_filename));*/
      for (file_idx = 0;
           file_idx < unit_info->num_of_files && !had_err && !found;
           file_idx++)
      {
        GtStr *encseq_filename;
        char *encseq_basename;
        GtStr *basename = NULL;

        encseq_filename = gt_str_array_get_str(unit_info->file_names, file_idx);
        encseq_basename = gt_basename(gt_str_get(encseq_filename));
        basename = gt_str_new_cstr(encseq_basename);
        gt_free(encseq_basename);
        if (0 == gt_str_cmp(basename,
                            mapping_filename))
        {
          if (file_set[file_idx] != false)
          {
            gt_error_set(err, "file %s double entry",
                         gt_str_get(mapping_filename));
            had_err = 1;
          }
          else
          {
            file_set[file_idx] = true;
            unit_info->map_files[file_idx] = genome_idx;
            found = 1;
            files_added++;
          }
        }
        gt_str_delete(basename);
      }
      if (!found && !had_err)
      {
        gt_error_set(err, "file %s not found in index, part of genome %s",
                     gt_str_get(mapping_filename),
                     gt_str_get(gt_str_array_get_str(unit_info->genome_names,
                                                     genome_idx)));
        had_err = 1;
      }
      found = 0;
      /*gt_lua_stack_dump(L);*/
      lua_pop(L, 1);
      gt_str_delete(mapping_filename);
    }
    lua_pop(L, 1);
    genome_idx++;
  }
  lua_pop(L, 1);

  gt_free(file_set);

  if (!had_err)
  {
    unit_info->num_of_genomes = genome_idx;
    /*fprintf(stderr, "num_of_genomes = %lu, genome_names size = %lu",
            unit_info->num_of_genomes,
            gt_str_array_size(unit_info->genome_names));*/
    gt_assert(unit_info->num_of_genomes ==
        gt_str_array_size(unit_info->genome_names));
  }
  if (!had_err && files_added != unit_info->num_of_files)
  {
    had_err = 1;
    gt_error_set(err, "different number of files in index"
                      " (%lu) than in unitfile (%lu)!",
                 unit_info->num_of_files,
                 files_added);
  }
  return had_err;
}

int gt_read_genomediff_unitfile(GtStr *unitfile,
                                struct GtShuUnitFileInfo_tag *unit_info,
                                GT_UNUSED GtLogger *logger,
                                GtError *err)
{
  int had_err = 0;
  /* open Lua */
  lua_State *L = luaL_newstate();
  /* library needed?
  lua_openlibs(L); */

  had_err = load(L, gt_str_get(unitfile), err);
  if (!had_err)
  {
    had_err = traverse_units(L, unit_info, err);
  }
  lua_close(L);
  return had_err;
}

void gt_delete_unit_file_info(struct GtShuUnitFileInfo_tag *unit_info) {
  gt_free(unit_info->map_files);
  gt_str_array_delete(unit_info->genome_names);
  gt_free(unit_info);
}
