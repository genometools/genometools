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

#ifndef GENOME_NODE_ITERATOR_LUA_H
#define GENOME_NODE_ITERATOR_LUA_H

#include "lua.h"

/* exports the GenomeNodeIterator class to Lua:

   -- Returns a new genome node iterator which performs a depth-first traversel
   -- of <genome_node> (including <genome_node> itself).
   function genome_node_iterator_new(genome_node)

   -- Returns a noew genome node iterator wich iterates over all direct children
   -- of <genome_node> (without <genome_node> itself).
   function genome_node_iterator_new_direct(genome_node)

   -- Returns the next genome node for <genome_node_iterator> or nil.
   function genome_node_iterator:next()
*/
int luaopen_genome_node_iterator(lua_State*);

#endif
