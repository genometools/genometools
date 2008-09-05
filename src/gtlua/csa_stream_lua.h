/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef CSA_STREAM_LUA_H
#define CSA_STREAM_LUA_H

#include "lua.h"

/* exports the CSAStream class (which implements the GenomeStream) interface
   to Lua:

   -- Returns a new CSA (consensus spliced alignment) stream object (a genome
   -- stream) which uses genome stream <in_stream> as input.
   -- The CSA stream replaces spliced alignments with computed consensus spliced
   -- alignments.
   -- The optional <join> parameters sets the length for the spliced alignment
   -- clustering (default: 300).
   function csa_stream_new(in_stream, join)
*/
int luaopen_csa_stream(lua_State*);

#endif
