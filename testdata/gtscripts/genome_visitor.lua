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

-- testing the Lua bindings for the GFF3 GenomeVisitor (similar to the gff3
-- tool)

function usage()
  io.stderr:write(string.format("Usage: %s GFF3_file\n", arg[0]))
  io.stderr:write("Parse and output the given GFF3_file.\n")
  os.exit(1)
end

if #arg == 1 then
  gff3file = arg[1]
else
  usage()
end

gs = gt.gff3_in_stream_new_sorted(gff3file)
gv = gt.gff3_visitor_new()
gn = gs:next_tree()
while (gn) do
  gn:accept(gv)
  gn = gs:next_tree()
end
