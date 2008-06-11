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

Testsuite = {}

function Testsuite:new(keywords)
  assert(keywords)
  o = {}
  o.keywords = keywords
  setmetatable(o, self)
  self.__index = self
  return o
end

function Testsuite:add_content(content)
  assert(content)
  self.contents = self.contents or {}
  self.contents[#self.contents + 1] = content
end

function Testsuite:add_test(name, def)
  assert(name and def)
  self.tests = self.tests or {}
  local test = {}
  test.name = name
  test.def = def
  self.tests[#self.tests + 1] = test
end

function Testsuite:add_directory(src, dest)
  assert(src and dest)
  self.directories = self.directories or {}
  local dir = {}
  dir.src = src
  dir.dest = dest
  self.directories[#self.directories + 1] = dir
end

function Testsuite:add_file(name, content)
  assert(name and content)
  self.files = self.files or {}
  local file = {}
  file.name = name
  file.content = content
  self.files[#self.files + 1] = file
end

function Testsuite:write(dir)
  assert(dir)
  if self.tests or self.contents then
    local outfile = io.open(dir .. "/testsuite.rb", "w")
    outfile:write([[
if $0 == __FILE__
  require 'stest'
  at_exit do
    OnError do exit 1 end
  end
end

$bin=File.join(Dir.pwd, "")
]])
    if self.tests then
      for _, test in ipairs(self.tests) do
        outfile:write("Name \"" .. test.name .. "\"\n")
        outfile:write("Keywords \"" .. self.keywords .. "\"\n")
        outfile:write("Test do\n");
        outfile:write(test.def)
        outfile:write("end\n\n")
      end
    end
    if self.contents then
      for _, content in ipairs(self.contents) do
        outfile:write("\n")
        outfile:write(content)
      end
    end
    outfile:close()
  end
  if self.directories then
    for _, directory in ipairs(self.directories) do
      assert(lfs.mkdir(dir .. "/" .. directory.dest))
      for file in lfs.dir(Extractor_project_home .. "/" .. directory.src) do
        if file ~= "." and file ~= ".." then
          print(directory.src .. "/" .. file)
          local infile = io.open(Extractor_project_home .. "/" ..
                                 directory.src .. "/" .. file)
          assert(infile, err)
          local outfile, err = io.open(dir .. "/" .. directory.dest .. "/" ..
                                       file, "w")
          assert(outfile, err)
          outfile:write(infile:read("*a"))
          infile:close()
          outfile:close()
        end
      end
    end
  end
  if self.files then
    for _, file in ipairs(self.files) do
      local outfile = io.open(dir .. "/" .. file.name, "w")
      outfile:write(file.content)
      outfile:close()
    end
  end
end
