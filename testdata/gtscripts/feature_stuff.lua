--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
]]

-- testing the Lua bindings for FeatureIndex and FeatureStream classes

function usage()
io.stderr:write(string.format("Usage: %s testdata_dir\n", arg[0]))
  io.stderr:write("Test the FeatureIndex and FeatureStream bindings.\n")
  os.exit(1)
end


if #arg == 1 then
  testdata = arg[1]
else
  usage()
end

-- set up the feature stream
genome_stream = gt.gff3_in_stream_new_sorted(testdata.."/gff3_file_1_short.txt")
feature_index = gt.feature_index_new()
genome_stream = gt.feature_stream_new(genome_stream, feature_index)
collectgarbage()

feature = genome_stream:next_tree()
while (feature) do
  feature = genome_stream:next_tree()
end

features = feature_index:get_features_for_seqid("ctg123");
assert(features)
gff3_visitor = gt.gff3_visitor()

for i,feature in ipairs(features) do
  feature:accept(gff3_visitor)
end
