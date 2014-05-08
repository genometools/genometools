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

#ifndef GENOME_NODE_LUA_H
#define GENOME_NODE_LUA_H

#include "lua.h"
#include "extended/genome_node.h"

/* exports the GenomeNode interface and its implementors to Lua:

   -- Returns the filename of <genome_node>.
   function genome_node:get_filename()

   -- Returns the range of <genome_node>.
   function genome_node:get_range()

   -- Returns the sequence id of <genome_node>.
   function genome_node:get_seqid()

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

*/

/* exports the FeatureNode class to Lua:

   -- Create a new feature node on sequence with ID <seqid> and type <type>
   -- which lies from <startpos> to <end> on strand <strand>.
   -- <startpos> and <endpos> always refer to the forward strand, therefore
   -- <startpos> has to bo smaller or equal than <endpos>.
   function feature_node_new(seqid, type, startpos, endpos, strand)

   -- Returns the strand of <feature_node>.
   function feature_node:get_strand()

   -- Returns the source of <feature_node>.
   function feature_node:get_source()

   -- Returns the score of <feature_node>.
   function feature_node:get_score()

   -- Returns the <attrib> attribute of <feature_node>.
   function feature_node:get_attribute(attrib)

   -- Returns an array containing the exons of <feature_node>.
   function feature_node:get_exons()

   -- Set the source of <feature_node> to <source>.
   function feature_node:set_source(source)

   -- Show leading part of GFF3 output for <feature_node>
   function feature_node:output_leading()

   -- Return type of <feature_node> as string.
   function feature_node:get_type()

   -- Extract the sequence of <feature_node>.
   -- If <join> is false and <feature_node> has type <type> the sequence is
   -- returned (using <region_mapping> to get it).
   -- If <join> is true and <feature_node> has children of type <type> their
   -- joined sequences are returned.
   -- If none of the above applies nil is returned.
   function feature_node:extract_sequence(type, join, region_mapping)

*/

/* exports the RegionNode class to Lua:

   -- Returns a new region node for sequence id <seqid> spanning <range>.
   function region_node_new(seqid, range)
*/

/* exports the MetaNode class to Lua:

   -- Returns a new region node with key <directive> and data string <data>.
   function meta_node_new(directive, data)

   -- Return directive of <meta_node> as string.
   function meta_node:get_directive()

   -- Return data of <meta_node> as string.
   function meta_node:get_data()
*/

int gt_lua_open_genome_node(lua_State*);

/* Push a <GtGenomeNode*> to Lua, takes ownership! */
void gt_lua_genome_node_push(lua_State*, GtGenomeNode*);

#define GENOME_NODE_METATABLE  "GenomeTools.genome_node"
#define check_genome_node(L, POS) \
                (GtGenomeNode**) luaL_checkudata(L, POS, GENOME_NODE_METATABLE)

#endif
