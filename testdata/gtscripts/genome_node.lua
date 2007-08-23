--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
]]

-- testing the Lua bindings for the GenomeNode interface

-- testing gt.genome_feature_new
rval, err = pcall(gt.genome_feature_new, nil, 1, 100, "+")
assert(not rval)
rval, err = pcall(gt.genome_feature_new, "not_defined", 1, 100, "+")
assert(not rval)
assert(string.find(err, "invalid feature type"))
rval, err = pcall(gt.genome_feature_new, "gene", "test", 100, "+")
assert(not rval)
rval, err = pcall(gt.genome_feature_new, "gene", 1, "test", "+")
assert(not rval)
rval, err = pcall(gt.genome_feature_new, "gene", 100, 1, "+")
assert(not rval)
assert(string.find(err, "must be <= end"))
rval, err = pcall(gt.genome_feature_new, "gene", 1, 100, "plus")
assert(not rval)
assert(string.find(err, "strand string must have length 1"))
rval, err = pcall(gt.genome_feature_new, "gene", 1, 100, "p")
assert(not rval)
assert(string.find(err, "invalid strand"))
gn = gt.genome_feature_new("gene", 1, 100, "+")

-- testing genome_node:get_filename
rval, fn = pcall(gn.get_filename, gn)
assert(rval)
assert(string.find(fn, "^Lua$"))
