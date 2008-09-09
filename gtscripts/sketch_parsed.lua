function usage()
  io.stderr:write(string.format("Usage: %s PNG_file GFF3_file\n", arg[0]))
  io.stderr:write("Create PNG representation of GFF3 annotation file.\n")
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
in_stream = nil; collectgarbage()
gn = feature_stream:next_tree()
-- fill feature index
while (gn) do
  gn = feature_stream:next_tree()
end

-- create diagram for first sequence ID in feature index
seqid = feature_index:get_first_seqid()
range = feature_index:get_range_for_seqid(seqid)
diagram = gt.diagram_new(feature_index, seqid, range)

-- create canvas
canvas = gt.canvas_new_png(600, nil)

-- sketch diagram on canvas
diagram:sketch(canvas)
-- write canvas to file
canvas:to_file(pngfile)
