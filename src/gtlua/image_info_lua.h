/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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

#ifndef IMAGE_INFO_LUA_H
#define IMAGE_INFO_LUA_H

#include "lua.h"

/* exports the Imageinfo class to Lua:

   -- returns a new ImageInfo object.
   function imageinfo_new()

   -- returns an array of tables with the fields "nw_x","nw_y","se_x","se_y" and
   -- "feature_ref" with the top left and bottom right coordinates in pixels or
   -- points and a GenomeNode reference per element drawn.
   function imageinfo:get_recmaps()
*/
int gt_lua_open_imageinfo(lua_State*);

#define IMAGEINFO_METATABLE  "GenomeTools.imageinfo"
#define check_imageinfo(L, POS) \
              (GtImageInfo**) luaL_checkudata(L, POS, IMAGEINFO_METATABLE)

#endif
