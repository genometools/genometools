--[[
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

-- testing the Lua bindings for the CSAStream binding

function usage()
  io.stderr:write(string.format("Usage: %s Fasta_file GFF3_file\n", arg[0]))
  io.stderr:write("Add CDS features to exon features given in GFF3_file.\n")
  os.exit(1)
end

if #arg == 2 then
  fasta_file = arg[1]
  gff3file = arg[2]
else
  usage()
end

in_stream = gt.gff3_in_stream_new_sorted(gff3file)
region_mapping = gt.region_mapping_new_seqfile(fasta_file)
cds_stream = gt.cds_stream_new(in_stream, region_mapping)
out_stream = gt.gff3_out_stream_new(cds_stream)
in_stream = nil;
cds_stream = nil;
collectgarbage() -- being nasty

gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end
