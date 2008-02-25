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

#ifndef GENOME_NODE_LUA_H
#define GENOME_NODE_LUA_H

#include "lua.h"
#include "libgtext/genome_node.h"

/* exports the GenomeNode interface and its implementors to Lua:

   -- Returns a new genome feature of <type> spanning <range> on <strand>.
   function genome_feature_new(type, range, strand)

   -- Returns a new sequence region for sequence id <seqid> spanning <range>.
   function sequence_region_new(seqid, range)

   -- Returns the filenname of <genome_node>.

   function genome_node:get_filename()

   -- Returns the range of <genome_node>.
   function genome_node:get_range()

   -- Returns the sequence id of <genome_node>.
   function genome_node:get_seqid()

   -- Set the sequence id of <genome_node> to <seqid>.
   function genome_node:set_seqid(seqid)

   -- Returns the strand of <genome_feature>.
   function genome_feature:get_strand()

   -- Returns the source of <genome_feature>.
   function genome_feature:get_source()

   -- Set the source of <genome_feature> to <source>.
   function genome_feature:set_source(source)

   -- Accept <genome_visitor>.
   function genome_node:accept(genome_visitor)

   -- Make <genome_node> the parent of <child_node>.
   function genome_node:is_part_of_genome_node(child_node)

   -- Mark <genome_node>.
   function genome_node:mark()

   -- Returns true if <genome_node> is marked, false otherwise.
   function genome_node:is_marked()

   -- Returns true if <genome_node> contains a marked node, false otherwise.
   function genome_node:contains_marked()

   -- Show leading part of GFF3 output for <genome_feature>
   function genome_feature:output_leading()

   -- Return type of <genome_feature> as string.
   function genome_feature:get_type()

   -- Extract the sequence of <genome_feature>.
   -- If <join> is false and <genome_feature> has type <type> the sequence is
   -- returned (using <region_mapping> to get it).
   -- If <join> is true and <genome_feature> has children of type <type> their
   -- joined sequences are returned.
   -- If none of the above applies nil is returned.
   function genome_feature:extract_sequence(type, join, region_mapping)
*/
int luaopen_genome_node(lua_State*);

/* push a GenomeNode to Lua, takes ownership! */
void genome_node_lua_push(lua_State*, GenomeNode*);

#define GENOME_NODE_METATABLE  "GenomeTools.genome_node"
#define check_genome_node(L, POS) \
                (GenomeNode**) luaL_checkudata(L, POS, GENOME_NODE_METATABLE)

#endif
