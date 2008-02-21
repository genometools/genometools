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

#include <assert.h>
#include "libgtlua/cds_stream_lua.h"
#include "libgtlua/csa_stream_lua.h"
#include "libgtlua/genome_node_lua.h"
#include "libgtlua/genome_node_iterator_lua.h"
#include "libgtlua/genome_stream_lua.h"
#include "libgtlua/genome_visitor_lua.h"
#include "libgtlua/gtext_lua.h"
#include "libgtlua/region_mapping_lua.h"
#include "libgtlua/stream_evaluator_lua.h"

int luaopen_gtext(lua_State *L)
{
  assert(L);
  luaopen_cds_stream(L);
  luaopen_csa_stream(L);
  luaopen_genome_node(L);
  luaopen_genome_node_iterator(L);
  luaopen_genome_stream(L);
  luaopen_genome_visitor(L);
  luaopen_region_mapping(L);
  luaopen_stream_evaluator(L);
  return 1;
}
