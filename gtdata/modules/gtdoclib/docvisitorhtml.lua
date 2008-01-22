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

lp = require 'cgilua/lp'

DocVisitorHTML = {}

local template_dir

function DocVisitorHTML:new(template_path)
  assert(template_path)
  template_dir = template_path
  o = {}
  setmetatable(o, self)
  self.__index = self
  return o
end

local function include(template, env)
  assert(template)
  local template_path = template_dir .. template
  env = env or {}
  env.io = io
  env.os = os
  env.ipairs = ipairs
  return lp.include(template_path, env)
end

local function codify(str)
  assert(str)
  local res = string.gsub(str, "<(.-)>", "<code>%1</code>")
  return string.gsub(res, " ([%a_][%a%d_%.]-%(%))", " <code>%1</code>")
end

function DocVisitorHTML:show_header()
  include("header.lp")
end

function DocVisitorHTML:visit_classes(classes)
  assert(classes)
  include("classes.lp", { classes = classes })
end

function DocVisitorHTML:visit_class(classname)
  assert(classname)
  include("class.lp", { classname = classname })
end

local sole_function_visited = false

function DocVisitorHTML:visit_sole_function(desc)
  if not sole_function_visited then
    include("sole_function.lp")
    sole_function_visited = true
  end
  self:visit_method(desc)
end

function DocVisitorHTML:visit_method(desc)
  assert(desc)
  include("method.lp", { name = desc.name, args = desc.args,
                         comment = codify(desc.comment) })
end

function DocVisitorHTML:show_footer()
  include("footer.lp")
end
