/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef STREAM_EVALUATOR_LUA_H
#define STREAM_EVALUATOR_LUA_H

#include "lua.h"

/* exports the StreamEvaluator class to Lua:

   stream_evaluator = gt.stream_evaluator_new(reality_stream, prediction_stream)
                      stream_evaluator:evaluate([genome_visitor])
                      stream_evaluaotr:show()
*/
int luaopen_stream_evaluator(lua_State*);

#endif
