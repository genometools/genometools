--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
]]

-- testing the Lua bindings for the GenomeStream interface

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

-- testing gt.gff3_in_stream_new_sorted
rval, err = pcall(gt.gff3_in_stream_new_sorted, "undefined")
assert(not rval)
assert(string.find(err, "does not exist"))

-- test correct file
gs = gt.gff3_in_stream_new_sorted(testdata.."/gff3_file_1_short.txt")
gn = gs:next_tree()
while (gn) do
  -- do something with the node...
  gn = gs:next_tree()
end

-- test corrupt file
gs = gt.gff3_in_stream_new_sorted(testdata.."/gt_gff3_fail_1.gff3")
gn = gs:next_tree()
while (gn) do
  rval, err = pcall(gs.next_tree, gs)
  if not rval then break end
end
assert(string.find(err, "already been defined"))

-- test unsorted file
gs = gt.gff3_in_stream_new_sorted(testdata.."/unsorted_gff3_file.txt")
gn = gs:next_tree()
while (gn) do
  rval, err = pcall(gs.next_tree, gs)
  if not rval then break end
end
assert(string.find(err, "is not sorted"))
