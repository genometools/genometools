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

#ifndef SCORE_MATRIX_LUA_H
#define SCORE_MATRIX_LUA_H

#include "lua.h"

/* exports the ScoreMatrix class to Lua:

   -- Returns a new protein score matrix object which has been read from file
   -- <path>.
   function score_matrix_new_read_protein(path)

   -- Returns the dimension of the <score_matrix> as number.
   function score_matrix:get_dimension()

   -- Returns the score for <idx1>,<idx2> as number.
   function score_matrix:get_score(idx1, idx2)
*/
int luaopen_score_matrix(lua_State*);

#endif
