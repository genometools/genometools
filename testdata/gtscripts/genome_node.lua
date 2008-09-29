--[[
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
]]

-- testing the Lua bindings for the GenomeNode interface

-- testing gt.feature_node_new
range = gt.range_new(1, 100)
rval, err = pcall(gt.feature_node_new, nil, nil, range:get_start(), range:get_end(), "+")
assert(not rval)
rval, err = pcall(gt.feature_node_new, "seqid", nil, range:get_start(), range:get_end(), "+")
assert(not rval)
rval, err = pcall(gt.feature_node_new, "seqid", "gene", "test", "+")
assert(not rval)
rval, err = pcall(gt.feature_node_new, "seqid", "gene", range:get_start(), range:get_end(), "plus")
assert(not rval)
assert(string.find(err, "strand string must have length 1"))
rval, err = pcall(gt.feature_node_new, "seqid", "gene", range:get_start(), range:get_end(), "p")
assert(not rval)
assert(string.find(err, "invalid strand"))
gn = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
assert(not gn:is_marked())
gn:mark()
assert(gn:is_marked())

parent = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
child  = gt.feature_node_new("seqid", "exon", range:get_start(), range:get_end(), "+")
parent:add_child(child)
assert(not parent:is_marked(parent))
assert(not parent:contains_marked(parent))
child:mark()
child  = nil; collectgarbage() -- being nasty
assert(not parent:is_marked(parent))
assert(parent:contains_marked(parent))

-- testing genome_node:get_filename
rval, fn = pcall(gn.get_filename, gn)
assert(rval)
assert(string.find(fn, "^generated$"))

-- testing gt.region_node_new
range = gt.range_new(1, 100)
rval, err = pcall(gt.region_node_new, nil, range:get_start(), range:get_end())
assert(not rval)
rval, err = pcall(gt.region_node_new, "chr1", "test")
assert(not rval)
gn = gt.region_node_new("chr1", range:get_start(), range:get_end())
