--[[
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

-- testing the Lua bindings for AnnotationSketch (similar to the sketch tool)

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
feature_index = gt.feature_index_memory_new()
feature_stream = gt.feature_stream_new(in_stream, feature_index)
in_stream = nil; collectgarbage() -- being nasty
gn = feature_stream:next_tree()
-- fill feature index
while (gn) do
  gn = feature_stream:next_tree()
end

seqid = feature_index:get_first_seqid()
range = feature_index:get_range_for_seqid(seqid)

ii = gt.imageinfo_new()

diagram = gt.diagram_new(feature_index, seqid, range)
layout = gt.layout_new(diagram, 800)
height = layout:get_height()
canvas = gt.canvas_cairo_file_new_png(800, height, ii)
-- nil as ImageInfo parameter must be ok
canvas2 = gt.canvas_cairo_file_new_png(800, height, nil)

-- test sketching with image info
layout:sketch(canvas)
canvas:to_file(pngfile)
if #(ii:get_recmaps()) ~= 16 then
  os.exit(1)
end

-- and without
layout:sketch(canvas2)
canvas2:to_file(pngfile)
