/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "extended/cds_stream_api.h"
#include "extended/luahelper.h"
#include "gtlua/cds_stream_lua.h"
#include "gtlua/genome_stream_lua.h"
#include "gtlua/region_mapping_lua.h"

static int cds_stream_lua_new(lua_State *L)
{
  GtNodeStream **cds_stream, **in_stream;
  GtRegionMapping **region_mapping;
  in_stream = check_genome_stream(L, 1);
  region_mapping = check_region_mapping(L, 2);
  cds_stream = lua_newuserdata(L, sizeof (GtNodeStream*));
  gt_assert(cds_stream);
  *cds_stream = gt_cds_stream_new(*in_stream,
                                  gt_region_mapping_ref(*region_mapping),
                                  1, /* XXX */
                                  "gtscript",
                                  true, /* XXX */
                                  true, /* XXX */
                                  false /* XXX */);
  luaL_getmetatable(L, GENOME_STREAM_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static const struct luaL_Reg cds_stream_lib_f [] = {
  { "cds_stream_new", cds_stream_lua_new },
  { NULL, NULL }
};

int gt_lua_open_cds_stream(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_register(L, "gt", cds_stream_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}
