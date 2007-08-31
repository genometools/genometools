--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
]]

-- testing the Lua bindings for libgtview (similar to the view tool)

function usage()
  io.stderr:write(string.format("Usage: %s PNG_file GFF3_file\n", arg[0]))
  io.stderr:write("Parse and output the given GFF3_file.\n")
  os.exit(1)
end

if #arg == 2 then
  pngfile  = arg[1]
  gff3file = arg[2]
else
  usage()
end

in_stream = gt.gff3_in_stream_new_sorted(gff3file)
feature_index = gt.feature_index_new()
feature_stream = gt.feature_stream_new(in_stream, feature_index)
in_stream = nil; collectgarbage() -- being nasty
gn = feature_stream:next_tree()
-- fill feature index
while (gn) do
  gn = feature_stream:next_tree()
end

seqid = feature_index:get_first_seqid()
startpos, endpos = feature_index:get_range_for_seqid(seqid)

diagram = gt.diagram_new(feature_index, startpos, endpos, seqid)
render = gt.render_new()
render:to_png(diagram, pngfile)
