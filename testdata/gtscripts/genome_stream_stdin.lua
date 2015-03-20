#!/usr/bin/env gt

instream = gt.gff3_in_stream_new_sorted()
out_stream = gt.gff3_out_stream_new(instream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end

instream = gt.gff3_in_stream_new_sorted(arg[1])
out_stream = gt.gff3_out_stream_new(instream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end