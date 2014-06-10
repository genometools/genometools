/*
  Copyright (c) 2008 Gordon Gremme <gordon@gremme.org>
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

#ifndef REGION_MAPPING_LUA_H
#define REGION_MAPPING_LUA_H

#include "lua.h"
#include "extended/region_mapping_api.h"

/* exports the RegionMapping class to Lua:

   -- Returns a new region mapping which maps everything onto sequence file
   -- <seqfile>.
   function region_mapping_new_seqfile(seqfile)
*/
int gt_lua_open_region_mapping(lua_State*);

/* Push a <GtRegionMapping*> to Lua, takes ownership! */
void gt_lua_region_mapping_push(lua_State *L, GtRegionMapping *rm);

#define REGION_MAPPING_METATABLE  "GenomeTools.region_mapping"
#define check_region_mapping(L, POS) \
        (GtRegionMapping**) luaL_checkudata(L, POS, REGION_MAPPING_METATABLE)

#endif
