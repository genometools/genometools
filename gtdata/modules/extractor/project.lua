--[[
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

require 'lfs'

Project = { name = "undefinded_name" }

function Project:new(home)
  assert(home)
  o = {}
  o.home = home
  -- store home in global environment
  _G.Extractor_project_home = home
  setmetatable(o, self)
  self.__index = self
  return o
end

function Project:set_name(name)
  assert(name)
  self.name = name
end

function Project:add(file)
  assert(file)
  self.parts = self.parts or {}
  self.parts[#self.parts + 1] = file
end

function Project:add_stest()
  local stest_file = File:new("testsuite/stest.rb")
  p:add(stest_file)
  stest_file = File:new("testsuite/stest_tests.rb")
  p:add(stest_file)
end

function Project:set_makefile(makefile)
  assert(makefile)
  self.makefile = makefile
end

function Project:write_tar_file()
  assert(lfs.mkdir(self.name))
  -- write makefile if necessary
  if self.makefile then
    self.makefile:write(self.name)
  end
  -- write other files
  for _, file in ipairs(self.parts) do
    file:write(self.name)
  end
  -- create archive
  local archivename = self.name .. ".tar.gz"
  print("writing archive '" .. archivename .. "'")
  os.execute("tar cvzf " .. archivename .. " " .. self.name)
end
