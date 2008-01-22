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
  o.functions = {}
  setmetatable(o, self)
  self.__index = self
  return o
end

function DocBase:add_class(classname, be_verbose)
  assert(classname)
  if be_verbose then
    print("class added: " .. classname)
  end
  self.classes[#self.classes + 1] = classname
end

function DocBase:add_function(funcname, funcargs, comment, be_verbose)
  assert(funcname and funcargs and comment)
  local desc = {}
  desc.name = funcname
  desc.args = funcargs
  desc.comment = comment
  if be_verbose then
    print("function added: " .. funcname)
  end
  self.functions[#self.functions + 1] = desc
end

local function function_keyword(ast)
  for i, keyword in ipairs(ast) do
    if keyword == "function" then
      return i
    end
  end
  return 0
end

function DocBase:process_ast(ast, be_verbose)
  assert(ast)
  for _, v in ipairs(ast) do
    if type(v) == "table" then
      self:process_ast(v, be_verbose)
    else
      local keyword = ast[1]
      -- print(keyword)
      if keyword == "class" then
        assert(#ast == 2)
        self["add_" .. ast[1]](self, ast[2], be_verbose)
        break
      elseif keyword == "comment" then
        local funcpos = function_keyword(ast)
        local complete_comment = ""
        if funcpos > 0 then
          assert(funcpos > 2)
          assert(#ast == funcpos + 2)
          if ast[2] == "undefined" then
            print("warning: undefined comment") -- XXX
          else
            complete_comment = table.concat(ast, " ", 2, funcpos-1)
          end
            self:add_function(ast[funcpos+1], ast[funcpos+2], complete_comment,
                              be_verbose)
          break
        end
      end
    end
  end
end

function DocBase:accept(visitor)
  assert(visitor)
  -- visit classes
  for _, classname in ipairs(self.classes) do
    visitor:visit_class(classname)
  end
  -- visit functions
  for _, funcdesc in ipairs(self.functions) do
    visitor:visit_function(funcdesc)
  end
end
