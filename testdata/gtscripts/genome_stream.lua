--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
]]

-- testing the Lua bindings for the GenomeStream interface

-- testing gt.gff3_in_stream_new_sorted
rval, err = pcall(gt.gff3_in_stream_new_sorted, "undefined")
assert(not rval)
assert(string.find(err, "does not exist"))
