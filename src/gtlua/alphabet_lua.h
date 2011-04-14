/*
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef ALPHABET_LUA_H
#define ALPHABET_LUA_H

#include "lua.h"
#include "core/alphabet.h"

/* exports the Alphabet class to Lua:

   -- Return a new protein alphabet.
   function alphabet_new_protein()

   -- Return a new DNA alphabet.
   function alphabet_new_dna()

   -- Return an empty alphabet.
   function alphabet_new_empty()

   -- Add the mapping of all given <characters> to the given <alphabet>.
   -- The first character is the result of subsequent <alphabet:decode()> calls.
   function alphabet:add_mapping(characters)

   -- Add <wildcard> to <alphabet>.
   function alphabet:add_wildcard(characters)

   -- Return a string containing the decoded character of the <code> number.
   function alphabet:decode(code)

   -- Return the size of <alphabet> as a number.
   function alphabet:size()
*/
int gt_lua_open_alphabet(lua_State*);

/* Push a <GtAlphabet*> to Lua, takes ownership! */
void gt_lua_alphabet_push(lua_State*, GtAlphabet*);

#define ALPHABET_METATABLE  "GenomeTools.alphabet"
#define check_alphabet(L, POS) \
          (GtAlphabet**) luaL_checkudata(L, POS, ALPHABET_METATABLE)

#endif
