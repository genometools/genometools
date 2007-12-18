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

require 'extractor'

if #arg == 1 then
  gt_home = arg[1]
else
  io.stderr:write(string.format("Usage: %s gt_home\n", arg[0]))
  io.stderr:write("Extract swalign program from GenomeTools home directory " ..
                  "gt_home.\n")
  os.exit(1)
end

p = Project:new(gt_home)
p:set_name("swalign")

prog = Program:new()
p:add(prog)

m = Module:new("src/libgtcore/xansi")
m:bare_includes()
p:add(m)

mf = Makefile:new()
p:set_makefile(mf)


p:write_tar_file()
