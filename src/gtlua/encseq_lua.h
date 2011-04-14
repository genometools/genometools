/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#ifndef ENCSEQ_LUA_H
#define ENCSEQ_LUA_H

#include "lua.h"
#include "core/encseq_api.h"

/* exports the Encseq class to Lua:      (XXX: document me)

   -- Returns a XXX with XXX.
   function ...(...)

*/
int gt_lua_open_encseq(lua_State*);

/* Push a <GtEncseq*> to Lua, takes ownership! */
void gt_lua_encseq_push(lua_State*, GtEncseq*);

#define ENCSEQ_METATABLE  "GenomeTools.encseq"
#define check_encseq(L, POS) \
        (GtEncseq**) luaL_checkudata(L, POS, ENCSEQ_METATABLE)

#define ENCSEQ_ENCODER_METATABLE  "GenomeTools.encseq_encoder"
#define check_encseq_encoder(L, POS) \
        (GtEncseqEncoder**) luaL_checkudata(L, POS, ENCSEQ_ENCODER_METATABLE)

#define ENCSEQ_LOADER_METATABLE  "GenomeTools.encseq_loader"
#define check_encseq_loader(L, POS) \
        (GtEncseqLoader**) luaL_checkudata(L, POS, ENCSEQ_LOADER_METATABLE)

#define ENCSEQ_BUILDER_METATABLE  "GenomeTools.encseq_builder"
#define check_encseq_builder(L, POS) \
        (GtEncseqBuilder**) luaL_checkudata(L, POS, ENCSEQ_BUILDER_METATABLE)

#define ENCSEQ_READER_METATABLE  "GenomeTools.encseq_reader"
#define check_encseq_reader(L, POS) \
        (GtEncseqReader**) luaL_checkudata(L, POS, ENCSEQ_READER_METATABLE)

#endif
