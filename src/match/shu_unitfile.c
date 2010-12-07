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
#include "core/encseq_api.h"
#include "core/fileutils_api.h"
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
                          GtStrArray *genome_names,
                          const GtEncseq *encseq,
                          GT_UNUSED GtError *err)
{
  int had_err = 0;
  const GtStrArray *filenames = gt_encseq_filenames(encseq);
  unsigned long file_idx = 0,
                file_check,
                num_of_files = gt_str_array_size(filenames);
  GtStr *encseq_filename,
        *mapping_filename;

  gt_assert(filenames);

  lua_pushnil(L); /*the first outer key*/
  while (lua_next(L, -2) != 0 && !had_err)
  {
    gt_str_array_add_cstr(genome_names, lua_tostring(L, -2));

    /* save file_idx to check for empty genomenames */
    file_check = file_idx;

    lua_pushnil(L); /* the first inner key */
    while (lua_next(L, -2) != 0 && !had_err)
    {
      if (num_of_files <= file_idx)
      {
        had_err = 1;
        gt_error_set(err, "more files in unitfile than in index!");
      }
      else
      {
        encseq_filename = gt_str_array_get_str(filenames, file_idx);
        mapping_filename = gt_str_new_cstr(lua_tostring(L, -1));
        if (gt_str_cmp(encseq_filename,
                        mapping_filename))
        {
          had_err = 1;
          gt_error_set(err, "files %s and %s do not match, check index!\n",
                       gt_str_get(encseq_filename),
                       gt_str_get(mapping_filename));
        }
        gt_lua_stack_dump(L);
        file_idx++;
      }
      lua_pop(L, 1);
    }
    lua_pop(L, 1);

    if (!had_err && file_check == file_idx)
    {
      had_err = 1;
      gt_error_set(err, "empty genome (no files): %s\n!",
                   lua_tostring(L, -2));
    }
  }
  lua_pop(L, 1);

  if (!had_err && file_idx != num_of_files)
  {
    had_err = 1;
    gt_error_set(err, "more files in index (%lu) than in unitfile(%lu)!",
                 num_of_files, file_idx + 1);
  }
  return had_err;
}

int gt_read_genomediff_unitfile(GtStr *unitfile,
                                const GtEncseq *encseq,
                                GtStrArray *genome_names,
                                GT_UNUSED unsigned long *num_of_genomes,
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
    had_err = traverse_units(L, genome_names, encseq, err);
  }
  return had_err;
}
