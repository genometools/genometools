--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

function usage()
  io.stderr:write(string.format("Usage: %s PNG_file reality_file " ..
                                "prediction_file\n", arg[0]))
  io.stderr:write("Evaluate prediction_file against reality_file and write " ..
                  "out PNG_file.\n")
  os.exit(1)
end

if #arg == 3 then
  png_file  = arg[1]
  reality_file = arg[2]
  prediction_file = arg[3]
else
  usage()
end

function render_to_png()
  local diagram = gt.diagram_new(feature_index, range, seqid)
  local render =  gt.render_new()
  render:to_png(diagram, png_file, width)
end

function show()
  os.execute("display "..png_file)
end

function render_and_show()
  render_to_png()
  show()
end

-- process input files
reality_stream = gt.gff3_in_stream_new_sorted(reality_file)
prediction_stream = gt.gff3_in_stream_new_sorted(prediction_file)
stream_evaluator = gt.stream_evaluator_new(reality_stream, prediction_stream)

feature_index = gt.feature_index_new()
feature_visitor = gt.feature_visitor_new(feature_index)

stream_evaluator:evaluate(feature_visitor)
stream_evaluator:show()

-- view results
seqid = feature_index:get_first_seqid()
range = feature_index:get_range_for_seqid(seqid)

diagram = gt.diagram_new(feature_index, range, seqid)
width = 800
render = gt.render_new()
render:to_png(diagram, png_file, width)
