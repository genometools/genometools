--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

-- testing gt.genome_feature_new
range = gt.range_new(1, 100)
rval, err = pcall(gt.genome_feature_new, nil, range, "+")
assert(not rval)
rval, err = pcall(gt.genome_feature_new, "not_defined", range, "+")
assert(not rval)
assert(string.find(err, "invalid feature type"))
rval, err = pcall(gt.genome_feature_new, "gene", "test", "+")
assert(not rval)
assert(string.find(err, "range expected"))
rval, err = pcall(gt.genome_feature_new, "gene", range, "plus")
assert(not rval)
assert(string.find(err, "strand string must have length 1"))
rval, err = pcall(gt.genome_feature_new, "gene", range, "p")
assert(not rval)
assert(string.find(err, "invalid strand"))
gn = gt.genome_feature_new("gene", range, "+")

-- testing genome_node:get_filename
rval, fn = pcall(gn.get_filename, gn)
assert(rval)
assert(string.find(fn, "^Lua$"))
