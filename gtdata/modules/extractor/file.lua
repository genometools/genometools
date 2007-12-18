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

module(..., package.seeall)

File = {}

function File:new(filename, do_not_read)
  o = {}
  if filename then
    o.filename = filename
    o.basename = filename:match("^.*/(.+)$") or filename
    if not do_not_read then
      local infile, err = io.open(filename, "r")
      assert(infile, err)
      o.filecontent = infile:read("*a")
      infile:close()
    end
  end
  setmetatable(o, self)
  self.__index = self
  return o
end

function File:bare_includes()
  print("match")
  self.filecontent = self.filecontent:gsub('(\n#include ").-/(.-"\n)', '%1%2')
end

function File:write(dir)
  assert(dir)
  local outfile = io.open(dir .. "/" .. self.basename, "w")
  outfile:write(self.filecontent)
  outfile:close()
end
