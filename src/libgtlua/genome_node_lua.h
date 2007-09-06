/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GENOME_NODE_LUA_H
#define GENOME_NODE_LUA_H

#include "lua.h"
#include "libgtext/genome_node.h"

/* exports the GenomeNode interface and its implementors to Lua:

   genome_node = gt.genome_feature_new(type, range, strand)
   string      = genome_node:get_filename()
   range       = genome_node:get_range()
                 genome_node:accept(genome_visitor)
                 parent_node:is_part_of_genome_node(child_node)
                 genome_node:mark()
   boolean     = genome_node:is_marked()
   boolean     = genome_node:contains_marked()
*/
int luaopen_genome_node(lua_State*);

/* push a GenomeNode to Lua, takes ownership! */
void genome_node_lua_push(lua_State*, GenomeNode*);

#define GENOME_NODE_METATABLE  "GenomeTools.genome_node"

#endif
