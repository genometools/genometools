--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
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

reality_stream = gt.gff3_in_stream_new_sorted(reality_file)
prediction_stream = gt.gff3_in_stream_new_sorted(prediction_file)
stream_evaluator = gt.stream_evaluator_new(reality_stream, prediction_stream)

feature_index = gt.feature_index_new()
feature_visitor = gt.feature_visitor_new(feature_index)

stream_evaluator:evaluate(feature_visitor)
stream_evaluator:show()

seqid = feature_index:get_first_seqid()
startpos, endpos = feature_index:get_range_for_seqid(seqid)
features = feature_index:get_features_for_range(seqid, startpos, endpos)

diagram = gt.diagram_new(features, startpos, endpos)
render = gt.render_new()
render:to_png(diagram, png_file)
