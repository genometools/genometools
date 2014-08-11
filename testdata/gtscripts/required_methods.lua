function usage()
  io.stderr:write(string.format("Usage: %s <GFF annotation>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
  os.exit(1)
end

f_stream = gt.custom_stream_new_unsorted()
f_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
function f_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    node:get_range()
    node:get_seqid()
    node:get_filename()
    node:get_line_number()
  end
  return node
end

local gn = f_stream:next_tree()
while (gn) do
  gn = f_stream:next_tree()
end