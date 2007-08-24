--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
]]

-- testing the Lua bindings for the GenomeStream interface

--[[
function usage()
  io.stderr:write(string.format("Usage: %s testdata_dir\n", arg[0]))
  io.stderr:write("Test the GenomeStream bindings.\n")
  os.exit(1)
end

if arg[1] then
  testdata = arg[1]
else
  usage()
end
]]

-- testing gt.gff3_in_stream_new_sorted
rval, err = pcall(gt.gff3_in_stream_new_sorted, "undefined")
assert(not rval)
assert(string.find(err, "does not exist"))

--[[
gs = gt.gff3_in_stream_new_sorted(testdata.."/gff3_file_1_short.txt")
gn = gs:next_tree()
]]
