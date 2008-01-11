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

local function generate_add_function(name)
  Program["add_" .. name] = function(self, arg)
    assert(arg)
    self[name .. "s"] = self[name .. "s"] or {}
    self[name .. "s"][#self[name .. "s"] + 1] = arg
  end
end

generate_add_function("include")
generate_add_function("typedef")
generate_add_function("function")

function Program:add_define(defname, defcontent)
  assert(defname and defcontent)
  self.defines = self.defines or {}
  self.defines[#self.defines + 1] = defname .. "  " .. defcontent
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
  if self.defines then
    for _, def in ipairs(self.defines) do
      outfile:write("#define " .. def .. "\n")
    end
    outfile:write("\n")
  end
  if self.typedefs then
    for _, typedef in ipairs(self.typedefs) do
      outfile:write(typedef)
      outfile:write("\n")
    end
    outfile:write("\n")
  end
  if self.functions then
    for _, func in ipairs(self.functions) do
      outfile:write(func)
      outfile:write("\n")
    end
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
