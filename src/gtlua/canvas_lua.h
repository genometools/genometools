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

#ifndef CANVAS_LUA_H
#define CANVAS_LUA_H

#include "lua.h"

/* exports the Canvas class to Lua:

   -- Return a Canvas object which acts as a PNG drawing surface of
   -- width <width> to be passed to rendering functions as a visitor.
   -- An <imageinfo> object is filled with coordinate information if given.
   -- If not needed, pass nil as <imageinfo>.
   function canvas_cairo_file_new_png(width, imageinfo)

   -- Return a Canvas object which acts as a PDF drawing surface of
   -- width <width> to be passed to rendering functions as a visitor.
   -- An <imageinfo> object is filled with coordinate information if given.
   function canvas_cairo_file_new_pdf(width, imageinfo)

   -- Return a Canvas object which acts as a PS drawing surface of
   -- width <width> to be passed to rendering functions as a visitor.
   -- An <imageinfo> object is filled with coordinate information if given.
   function canvas_cairo_file_new_ps(width, imageinfo)

   -- Return a Canvas object which acts as a SVG drawing surface of
   -- width <width> to be passed to rendering functions as a visitor.
   -- An <imageinfo> object is filled with coordinate information if given.
   function canvas_cairo_file_new_svg(width, imageinfo)

   -- Creates an image file with the given <filename> which contains the
   -- contents of the canvas (only for objects created with
   -- <canvas_cairo_file_new_*()>).
   function canvas:to_file(filename)
*/
int gt_lua_open_canvas(lua_State*);

#define CANVAS_METATABLE  "GenomeTools.canvas"
#define check_canvas(L, POS) \
              (GtCanvas**) luaL_checkudata(L, POS, CANVAS_METATABLE)

#endif
