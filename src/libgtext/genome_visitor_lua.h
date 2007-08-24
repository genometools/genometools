/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_VISITOR_LUA_H
#define GENOME_VISITOR_LUA_H

#include "lua.h"

/* exports the GenomeVisitor interface and its implementors to Lua:

   genome_visitor = gt.gff3_visitor()
*/
int luaopen_genome_visitor(lua_State*);

#define GENOME_VISITOR_METATABLE  "GenomeTools.genome_visitor"
#define check_genome_visitor(L, POS) \
          (GenomeVisitor**) luaL_checkudata(L, POS, GENOME_VISITOR_METATABLE);

#endif
