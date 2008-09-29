require 'gtruby'

if ARGV.size != 3 then
  STDERR.puts "Usage: #{$0} style_file PNG_file GFF3_file"
  STDERR.puts "Create PNG representation of GFF3 annotation file."
  exit(1)
end

(stylefile, pngfile, gff3file) = ARGV

# load style file
style = GT::Style.new()
style.load_file(stylefile)

# create feature index
feature_index = GT::FeatureIndexMemory.new()

# add GFF3 file to index
feature_index.add_gff3file(gff3file)

# create diagram for first sequence ID in feature index
seqid = feature_index.get_first_seqid()
range = feature_index.get_range_for_seqid(seqid)
diagram = GT::Diagram.new(feature_index, seqid, range, style)

# create canvas
canvas = GT::CanvasCairoFile.new(style, 600, nil)

# sketch diagram on canvas
diagram.sketch(canvas)

# write canvas to file
canvas.to_file(pngfile)
