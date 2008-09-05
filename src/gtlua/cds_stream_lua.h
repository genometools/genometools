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

#ifndef CDS_STREAM_LUA_H
#define CDS_STREAM_LUA_H

#include "lua.h"

/* exports the CDSStream class (which implements the GenomeStream) interface
   to Lua:

   -- Returns a new CDS (coding sequence) stream object (a genome stream) which
   -- uses genome stream <in_stream> as input.
   -- The CDS stream adds CDS features to exon features in <in_stream>.
   -- The given <region_mapping> is used to map the sequence regions given in
   -- <in_stream> to the actual sequence files necessary for computing the
   -- coding sequences.
   function cds_stream_new(region_mapping)
*/
int luaopen_cds_stream(lua_State*);

#endif
