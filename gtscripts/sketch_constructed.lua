function usage()
  io.stderr:write(string.format("Usage: %s Style_file PNG_file\n", arg[0]))
  os.exit(1)
end

if #arg == 2 then
  style_file = arg[1]
  png_file   = arg[2]
else
  usage()
end

-- load style file
dofile(style_file)

-- construct the example features
seqid = "chromosome_21"
nodes = {}

-- construct a gene on the forward strand with two exons
gene   = gt.feature_node_new(seqid, "gene", 100, 900, "+")
exon   = gt.feature_node_new(seqid, "exon", 100, 200, "+")
gene:add_child(exon)
intron = gt.feature_node_new(seqid, "intron", 201, 799, "+")
gene:add_child(intron)
exon   = gt.feature_node_new(seqid, "exon", 800, 900, "+")
gene:add_child(exon)
nodes[1] = gene

-- construct a single-exon gene on the reverse strand
-- (within the intron of the forward strand gene)
reverse_gene = gt.feature_node_new(seqid, "gene", 400, 600, "-")
reverse_exon = gt.feature_node_new(seqid, "exon", 400, 600, "-")
reverse_gene:add_child(reverse_exon)
nodes[2] = reverse_gene

-- create diagram
diagram = gt.diagram_new_from_array(nodes, 1, 1000)
layout = gt.layout_new(diagram, 600)
height = layout:get_height()

-- create canvas
canvas = gt.canvas_cairo_file_new_png(600, height, nil)

-- sketch layout on canvas
layout:sketch(canvas)

-- write canvas to file
canvas:to_file(png_file)
