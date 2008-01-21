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

module(..., package.seeall)

DocBase = {}

function DocBase:new()
  o = {}
  o.classes = {}
  setmetatable(o, self)
  self.__index = self
  return o
end

function DocBase:add_class(classname)
  assert(classname)
  self.classes[#self.classes + 1] = classname
end

function DocBase:process_ast(ast)
  assert(ast)
  for _, v in ipairs(ast) do
    if type(v) == "table" then
      self:process_ast(v)
    else
      assert(#ast == 2)
      self["add_" .. ast[1]](self, ast[2])
      break
    end
  end
end

function DocBase:accept(visitor)
  assert(visitor)
  -- visit classes
  for _, classname in ipairs(self.classes) do
    visitor:visit_class(classname)
  end
end
