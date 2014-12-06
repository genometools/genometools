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

#ifndef FEATURE_NODE_ITERATOR_LUA_H
#define FEATURE_NODE_ITERATOR_LUA_H

#include "lua.h"

/* exports the FeatureNodeIterator class to Lua:

   -- Returns a new feature node iterator which performs a depth-first traversal
   -- of <node> (including <node> itself).
   function feature_node_iterator_new(node)

   -- Returns a new feature node iterator wich iterates over all direct children
   -- of <node> (without <node> itself).
   function feature_node_iterator_new_direct(node)

   -- Returns the next node for <feature_node_iterator> or nil.
   function feature_node_iterator:next()
*/
int gt_lua_open_feature_node_iterator(lua_State*);

#endif
