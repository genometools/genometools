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

require 'stringext'
local w = require 'warning'

DocBase = {}

function DocBase:new()
  o = {}
  o.classes   = {}
  o.solefuncs = {}
  setmetatable(o, self)
  self.__index = self
  return o
end

function DocBase:add_class(classname, be_verbose)
  assert(classname)
  if be_verbose then
    print("class added: " .. classname)
  end
  self.classes[classname] = self.classes[classname] or {}
end

function DocBase:add_method(funcname, funcargs, comment, be_verbose)
  assert(funcname and funcargs and comment)
  local desc = {}
  -- remove ``GenomeTools_'' prefix which is used to extend exported C classes
  desc.name = string.gsub(funcname, "^GenomeTools_", "")
  desc.args = funcargs
  desc.comment = comment
  if be_verbose then
    print("method added: " .. desc.name)
  end
  local classname
  -- check if function is a constructor
  classname = string.match(desc.name, "^(%a[%a%d_]*)_new[%a%d_]*")
  if not classname then
    -- check if function name starts with an upper-case letter and is a method
    classname = string.match(desc.name, "^(%a[%a%d_]*):")
  end
  -- transform classname, if necessary
  if classname then
    -- special case for abbrevated class names
    classname = string.gsub(classname, "^%a%a%a_", string.upper)
    classname = string.gsub(classname, "^%a", string.upper)
    classname = string.gsub(classname, "_%a", string.upper)
    classname = string.gsub(classname, "_", "")
  end
  if be_verbose and classname then
    print("classname found: " .. classname)
  end
  -- if this is a valid classname, try to store method in class
  if classname and self.classes[classname] then
    self.classes[classname][#self.classes[classname] + 1] = desc
  else
    self.solefuncs[#self.solefuncs + 1] = desc
  end
end

local function method_keyword(ast)
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
      if be_verbose then
        print("keyword: " .. keyword)
      end
      if keyword == "class" then
        assert(#ast == 2)
        self["add_" .. ast[1]](self, ast[2], be_verbose)
        break
      elseif keyword == "comment" then
        local funcpos = method_keyword(ast)
        local complete_comment = ""
        if funcpos > 0 then
          assert(funcpos > 2)
          assert(#ast == funcpos + 2)
          if ast[2] == "undefined" then
            w.warning("undefined comment")
          else
            complete_comment = table.concat(ast, "", 2, funcpos-1)
            complete_comment = string.strip(complete_comment)
          end
            self:add_method(ast[funcpos+1], ast[funcpos+2], complete_comment,
                            be_verbose)
          break
        end
      end
    end
  end
end

function DocBase:accept(visitor)
  assert(visitor)
  -- visit all classes
  local sorted_classes = {}
  for classname in pairs(self.classes) do
    if #self.classes[classname] > 0 then
      sorted_classes[#sorted_classes + 1] = classname
    end
  end
  table.sort(sorted_classes)
  if visitor.visit_classes then
    visitor:visit_classes(sorted_classes)
  end
  -- visit sole functions
  for _, funcdesc in ipairs(self.solefuncs) do
    if visitor.visit_sole_function then
      visitor:visit_sole_function(funcdesc)
    else
      visitor:visit_method(funcdesc)
    end
  end
  -- visit each class
  for _, classname in ipairs(sorted_classes) do
    visitor:visit_class(classname)
    -- visit methods for class
    for _, method in ipairs(self.classes[classname]) do
      visitor:visit_method(method)
    end
  end
end
