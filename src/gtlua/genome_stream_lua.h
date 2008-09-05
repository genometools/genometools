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

#ifndef GENOME_STREAM_LUA_H
#define GENOME_STREAM_LUA_H

#include "lua.h"

/* exports the GenomeStream interface and its implementors of libgtext to Lua:

   -- Returns a new GFF3 input stream object for <filename>. The file <filename>
   -- has to be a sorted GFF3 file.
   function gff3_in_stream_new_sorted(filename)

   -- Returns a new GFF3 output stream which pulls its features from
   -- <genome_stream>.
   function gff3_out_stream_new(genome_stream)

   -- Returns the next genome node for <genome_stream> or nil.
   function genome_stream:next_tree()
*/
int luaopen_genome_stream(lua_State*);

#define GENOME_STREAM_METATABLE  "GenomeTools.genome_stream"
#define check_genome_stream(L, POS) \
          (GenomeStream**) luaL_checkudata(L, POS, GENOME_STREAM_METATABLE)

#endif
