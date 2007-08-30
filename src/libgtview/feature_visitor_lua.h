/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FEATURE_VISITOR_LUA_H
#define FEATURE_VISITOR_LUA_H

#include "lua.h"

/* exports the FeatureVisitor class (which implements the GenomeVisitor)
   interface to Lua:

   genome_visitor = gt.feature_visitor_new(feature_index)
*/
int luaopen_feature_visitor(lua_State*);

#endif
