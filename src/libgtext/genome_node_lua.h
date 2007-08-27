/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_NODE_LUA_H
#define GENOME_NODE_LUA_H

#include "lua.h"
#include "libgtext/genome_node.h"

/* exports the GenomeNode interface and its implementors to Lua:

   genome_node = gt.genome_feature_new(type, startpos, endpos, strand)
   string      = genome_node:get_filename()
                 genome_node:accept(genome_visitor)
*/
int luaopen_genome_node(lua_State*);

/* push a GenomeNode to Lua, takes ownership! */
void genome_node_lua_push(lua_State*, GenomeNode*);

#define GENOME_NODE_METATABLE  "GenomeTools.genome_node"

#endif
