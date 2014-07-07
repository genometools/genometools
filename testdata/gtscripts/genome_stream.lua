--[[
  Copyright (c) 2007 Gordon Gremme <gordon@gremme.org>
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
rval, err = pcall(gs.next_tree, gs)
assert(not rval)
assert(string.find(err, "already been defined"))

-- test unsorted file
gs = gt.gff3_in_stream_new_sorted(testdata.."/unsorted_gff3_file.txt")
rval, err = pcall(gs.next_tree, gs)
assert(not rval)
assert(string.find(err, "is not sorted"))


-- test custom streams, missing override
cs = gt.custom_stream_new_sorted()
rval, err = pcall(cs.next_tree, cs)
assert(not rval)
assert(string.find(err, "method defined in custom stream"))

-- test custom streams, wrong return type
gs = gt.gff3_in_stream_new_sorted(testdata.."/eden.gff3")
cs = gt.custom_stream_new_sorted()
function cs:next_tree()
  if gs:next_tree() then
    return 1
  end
end
rval, err = pcall(gs.next_tree, cs)
assert(not rval)
assert(string.find(err, "return a genome node"))

-- test custom streams, runtime error
gs = gt.gff3_in_stream_new_sorted(testdata.."/eden.gff3")
cs = gt.custom_stream_new_sorted()
function cs:next_tree()
  n = gs:next_tree()
  x = 1 + nil
  return n
end
rval, err = pcall(gs.next_tree, cs)
assert(not rval)
assert(string.find(err, "perform arithmetic on a nil"))

-- test custom streams
gs = gt.gff3_in_stream_new_sorted(testdata.."/eden.gff3")
cs = gt.custom_stream_new_sorted()
cs.instream = gs
function cs:next_tree()
  return cs.instream:next_tree()
end
rval, err = pcall(gs.next_tree, cs)
assert(rval)

-- test custom streams
cs = gt.custom_stream_new_sorted()
cs.count = 1
function cs:next_tree()
  if cs.count <= 5 then
    n = gt.feature_node_new("test", "gene", 1+cs.count, 100+cs.count, "+")
    cs.count = cs.count + 1
  else
    n = nil
  end
  return n
end
nodes = {}
local gn = cs:next_tree()
while (gn) do
  table.insert(nodes, gn)
  gn = cs:next_tree()
end
assert(#nodes == 5)
