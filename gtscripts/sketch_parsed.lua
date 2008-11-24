function usage()
  io.stderr:write(string.format("Usage: %s Style_file PNG_file GFF3_file\n", arg[0]))
  io.stderr:write("Create PNG representation of GFF3 annotation file.\n")
  os.exit(1)
end

if #arg == 3 then
  style_file = arg[1]
  png_file   = arg[2]
  gff3_file  = arg[3]
else
  usage()
end

-- load style file
dofile(style_file)

-- create feature index
feature_index = gt.feature_index_memory_new()

-- add GFF3 file to index
feature_index:add_gff3file(gff3_file)

-- create diagram for first sequence ID in feature index
seqid = feature_index:get_first_seqid()
range = feature_index:get_range_for_seqid(seqid)
diagram = gt.diagram_new(feature_index, seqid, range)

-- create layout
layout = gt.layout_new(diagram, 600)
height = layout:get_height()

-- create canvas
canvas = gt.canvas_cairo_file_new_png(600, height, nil)

-- sketch layout on canvas
layout:sketch(canvas)

-- write canvas to file
canvas:to_file(png_file)
