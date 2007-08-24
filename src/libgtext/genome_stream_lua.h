/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_STREAM_LUA_H
#define GENOME_STREAM_LUA_H

#include "lua.h"

/* exports the GenomeStream interface and its implementors to Lua:

   genome_stream = gt.gff3_in_stream_new_sorted(filename)
   genome_stream = gt.gff3_out_stream_new(genome_stream)
   genome_node   = genome_stream:next_tree()
*/
int luaopen_genome_stream(lua_State*);

#endif
