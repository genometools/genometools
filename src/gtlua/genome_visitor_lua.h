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

#ifndef GENOME_VISITOR_LUA_H
#define GENOME_VISITOR_LUA_H

#include "lua.h"

/* exports the GtNodeVisitor interface and its implementors to Lua:

*/

/* exports the GFF3Visitor class to Lua:
   -- Returns a new GFF3 visitor.
   function gff3_visitor_new()
*/

int gt_lua_open_genome_visitor(lua_State*);

#define GENOME_VISITOR_METATABLE  "GenomeTools.genome_visitor"
#define check_genome_visitor(L, POS) \
             (GtNodeVisitor**) luaL_checkudata(L, POS, GENOME_VISITOR_METATABLE)

#endif
