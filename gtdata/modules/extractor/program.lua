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

Program = {}

function Program:new(progname)
  o = {}
  if progname then
    o.progname = progname .. ".c"
  else
    o.progname = "prog.c"
  end
  setmetatable(o, self)
  self.__index = self
  return o
end

function Program:add_include(inc)
  assert(inc)
  self.includes = self.includes or {}
  self.includes[#self.includes + 1] = inc
end

function Program:add_function(func)
  assert(func)
  self.functions = self.functions or {}
  self.functions[#self.functions + 1] = func
end

function Program:set_content(content)
  assert(content)
  self.content = content
end

function Program:write(dir)
  assert(dir)
  local outfile = io.open(dir .. "/" .. self.progname, "w")
  if self.includes then
    for _, inc in ipairs(self.includes) do
      outfile:write("#include " .. inc .. "\n")
    end
    outfile:write("\n")
  end
  for _, func in ipairs(self.functions) do
    outfile:write(func)
    outfile:write("\n")
  end
  outfile:write("int main(int argc, char *argv[])\n")
  outfile:write("{\n")
  if self.content then
    outfile:write(self.content)
  end
  outfile:write("  return 0;\n")
  outfile:write("}\n")
  outfile:close()
end
