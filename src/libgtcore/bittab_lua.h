/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef BITTAB_LUA_H
#define BITTAB_LUA_H

#include "lua.h"

/* exports the Bittab class to Lua:

   bittab = gt.bittab_new(num_of_bits)
   bittab:set_bit(bit)
   bittab:unset_bit(bit)
*/
int luaopen_bittab(lua_State*);

#endif
