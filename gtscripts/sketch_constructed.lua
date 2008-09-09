function usage()
  io.stderr:write(string.format("Usage: %s Style_file PNG_file\n", arg[0]))
  os.exit(1)
end

if #arg == 2 then
  stylefile = arg[1]
  pngfile   = arg[2]
else
  usage()
end

-- load style file
dofile(stylefile)

-- construct the example features
seqid = "chromosome_21"
nodes = {}

-- construct a gene on the forward strand with two exons
gene   = gt.genome_feature_new(seqid, "gene", gt.range_new(100, 900), "+")
exon   = gt.genome_feature_new(seqid, "exon", gt.range_new(100, 200), "+")
gene:is_part_of_genome_node(exon)
intron = gt.genome_feature_new(seqid, "intron", gt.range_new(201, 799), "+")
gene:is_part_of_genome_node(intron)
exon   = gt.genome_feature_new(seqid, "exon", gt.range_new(800, 900), "+")
gene:is_part_of_genome_node(exon)
nodes[1] = gene

-- construct a single-exon gene on the reverse strand
-- (within the intron of the forward strand gene)
reverse_gene = gt.genome_feature_new(seqid, "gene", gt.range_new(400, 600), "-")
reverse_exon = gt.genome_feature_new(seqid, "exon", gt.range_new(400, 600), "-")
reverse_gene:is_part_of_genome_node(reverse_exon)
nodes[2] = reverse_gene

-- create diagram
diagram = gt.diagram_new_from_array(nodes, gt.range_new(1, 1000))

-- create canvas
canvas = gt.canvas_new_png(600, nil)

-- sketch diagram on canvas
diagram:sketch(canvas)

-- write canvas to file
canvas:to_file(pngfile)
