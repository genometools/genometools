/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef BITTAB_LUA_H
#define BITTAB_LUA_H

#include "lua.h"

/* exports the Bittab class to Lua:

   -- Returns a bittab with <num_of_bits> many bits.
   function bittab_new(num_of_bits)

   -- Set <bit> in <bittab>.
   function bittab:set_bit(bit)

   -- Unset <bit> in <bittab>.
   function bittab:unset_bit(bit)

   -- Store the complement of bittab <src> in <bittab>.
   -- <bittab> and <src> must have the same size.
   function bittab:complement(src)

   -- Set <bittab> equal to bittab <src>.
   -- <bittab> and <src> must have the same size.
   function bittab:equal(src)

   -- Set <bittab> equal to the bitwise AND of <bittab> and <src>.
   -- <bittab> and <src> must have the same size.
   function bittab:and_equal(src)

   -- Returns true if <bit> is set in <bittab>, false otherwise.
   function bittab:bit_is_set(bit)
*/
int luaopen_bittab(lua_State*);

#endif
