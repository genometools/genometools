/*
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef RANGE_LUA_H
#define RANGE_LUA_H

#include "lua.h"
#include "core/range.h"

/* exports the Range class to Lua:

   -- Returns a new range object with start <startpos> and end <endpos>.
   -- <startpos> must be smaller or equal than <endpos>.
   function range_new(startpos, endpos)

   -- Returns start of <range>.
   function range:get_start()

   -- Returns end of <range>.
   function range:get_end()

   -- Returns true if <range> and <other_range> overlap, false otherwise.
   function range:overlap(other_range)

   -- Returns true if <range> and <other_range> overlap with delta <delta>,
   -- false otherwise.
   function range:overlap_delta(other_range, delta)

   -- Returns the length of <range>.
   function range:length()

   -- Returns true if <range> contains <other_range>, false otherwise.
   function range:contains(other_range)

   -- Returns a new range consisting of <range> and <other_range> joined.
   function range:join(other_range)

   -- Returns true if <point> lies within <range>, false otherwise.
   function range:within(point)

   -- Returns an array containing the ranges from array <range_array> in sorted
   -- order.
   function ranges_sort(range_array)

   -- Returns true if the ranges in array <range_array> are sorted, false
   -- otherwise.
   function ranges_are_sorted(range_array)
*/
int gt_lua_open_range(lua_State*);

/* push a GtRange to Lua, returns 1 (number of elements pushed) */
int gt_lua_range_push(lua_State*, GtRange);

#define RANGE_METATABLE  "GenomeTools.range"
#define check_range(L, POS) \
          (GtRange*) luaL_checkudata(L, POS, RANGE_METATABLE)

#endif
