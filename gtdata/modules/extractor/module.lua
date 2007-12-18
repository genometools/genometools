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

require "extractor.file"

Module = {}

function Module:new(modulename)
  assert(modulename)
  o = {}
  o.c_file = extractor.File:new(modulename .. ".c")
  o.h_file = extractor.File:new(modulename .. ".h")
  setmetatable(o, self)
  self.__index = self
  return o
end

function Module:bare_includes()
  self.c_file:bare_includes()
  self.h_file:bare_includes()
end

function Module:remove_include(header)
  self.c_file:remove_include(header)
  self.h_file:remove_include(header)
end

function Module:remove_function(func)
  assert(func)
  self.c_file:remove_function(func)
  self.h_file:remove_function(func)
end

function Module:remove_example()
  self.c_file:remove_example()
  self.h_file:remove_example()
end

function Module:remove_unit_test()
  self.c_file:remove_unit_test()
  self.h_file:remove_unit_test()
end

function Module:ma2xansi()
  self.c_file:ma2xansi()
  self.h_file:ma2xansi()
end

function Module:write(dir)
  assert(dir)
  self.c_file:write(dir)
  self.h_file:write(dir)
end
