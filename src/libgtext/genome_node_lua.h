/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_NODE_LUA_H
#define GENOME_NODE_LUA_H

#include "lua.h"

/* exports the GenomeNode interface and its implementors to Lua:

   genome_node = gt.genome_feature_new(type, start, end, strand)
   string      = genome_node:get_filename()
*/
int luaopen_genome_node(lua_State*);

#endif
