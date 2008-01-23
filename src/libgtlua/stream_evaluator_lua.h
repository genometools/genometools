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

#ifndef STREAM_EVALUATOR_LUA_H
#define STREAM_EVALUATOR_LUA_H

#include "lua.h"

/* exports the StreamEvaluator class to Lua:

   -- Returns a new stream evaluator object for the two genome streams
   -- <reality_stream> and <prediction_stream>.
   function stream_evaluator_new(reality_stream, prediction_stream)

   -- Run evaluation of <stream_evaluator>. All evaluated features are visited
   -- by the optional <genome_visitor>.
   function stream_evaluator:evaluate(genome_visitor)

   -- Show result of <stream_evaluator> on stdout.
   function stream_evaluator:show()
*/
int luaopen_stream_evaluator(lua_State*);

#endif
