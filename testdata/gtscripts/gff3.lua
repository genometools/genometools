--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
]]

-- testing the Lua bindings for the GFF3 output stream (similar to the gff3
-- tool)

function usage()
  io.stderr:write(string.format("Usage: %s GFF3_file\n", arg[0]))
  io.stderr:write("Parse and output the given GFF3_file.\n")
  os.exit(1)
end

if #arg == 1 then
  gff3file = arg[1]
else
  usage()
end

in_stream = gt.gff3_in_stream_new_sorted(gff3file)
out_stream = gt.gff3_out_stream_new(in_stream)
in_stream = nil; collectgarbage() -- being nasty
gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end
