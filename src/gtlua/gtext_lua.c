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

#include "core/assert_api.h"
#include "gtlua/cds_stream_lua.h"
#include "gtlua/csa_stream_lua.h"
#include "gtlua/feature_index_lua.h"
#include "gtlua/feature_node_iterator_lua.h"
#include "gtlua/feature_stream_lua.h"
#include "gtlua/feature_visitor_lua.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/genome_stream_lua.h"
#include "gtlua/genome_visitor_lua.h"
#include "gtlua/gtext_lua.h"
#include "gtlua/region_mapping_lua.h"
#include "gtlua/stream_evaluator_lua.h"

int gt_lua_open_extended(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  gt_lua_open_cds_stream(L);
  gt_lua_open_csa_stream(L);
  gt_lua_open_feature_index(L);
  gt_lua_open_feature_node_iterator(L);
  gt_lua_open_feature_stream(L);
  gt_lua_open_feature_visitor(L);
  gt_lua_open_genome_node(L);
  gt_lua_open_genome_stream(L);
  gt_lua_open_genome_visitor(L);
  gt_lua_open_region_mapping(L);
  gt_lua_open_stream_evaluator(L);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}
