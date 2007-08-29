/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FEATURE_STREAM_LUA_H
#define FEATURE_STREAM_LUA_H

#include "lua.h"

/* exports the FeatureStream class (which implements the GenomeStream) interface
   to Lua:

   feature_stream = gt.feature_stream_new(feature_index)
*/
int luaopen_feature_stream(lua_State*);

#endif
