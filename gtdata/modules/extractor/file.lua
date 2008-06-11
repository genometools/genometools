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

File = {}

function File:new(filename, do_not_read)
  o = {}
  if filename then
    o.filename = filename
    o.basename = filename:match("^.*/(.+)$") or filename
    if not do_not_read then
      local infile, err
      if _G.Extractor_project_home then
        infile, err = io.open(_G.Extractor_project_home .. "/" .. filename, "r")
      else
        infile, err = io.open(filename, "r")
      end
      assert(infile, err)
      o.filecontent = infile:read("*a")
      infile:close()
    end
  end
  setmetatable(o, self)
  self.__index = self
  return o
end

function File:get_typedef(typedef)
  assert(typedef)
  return self.filecontent:match("typedef struct .-" .. typedef .. ";")
end

function File:get_function(func)
  assert(func)
  return(self.filecontent:match('\n(%S*%s?%S+%s+' .. func ..
                                '%(.-%)\n{.-\n})'))
end

function File:bare_includes()
  self.filecontent = self.filecontent:gsub('(\n#include ").-/(.-")', '%1%2')
end

function File:remove_include(header)
  assert(header)
  self.filecontent = self.filecontent:gsub('#include "' .. header .. '"\n', '')
end

function File:remove_function(func)
  assert(func)
  if self.basename:match("%.h$") then -- header file
    -- somewhat hackish way to remove function name with preceding comment,
    -- because Lua doesn not support grouping patterns (to make the comment
    -- optional). If the comment is not optional the non-greedy comment matcher
    -- woudl ``eat to much''.
    replacement = {}
    replacement[func] = ''
    self.filecontent = self.filecontent:gsub('\n/%*.-%*/' ..
                                             '\n%S*%s?%S+%s+' .. '(%S+)' ..
                                             '%(.-%);', replacement)
    -- remove function name without preceding comment
    self.filecontent = self.filecontent:gsub('\n%S*%s?%S+%s+' .. func ..
                                             '%(.-%);', '')
  else -- C file
    self.filecontent = self.filecontent:gsub('\n%S*%s?%S+%s?%S+%s+' .. func ..
                                             '%(.-%)\n{.-\n}', '')
  end
  self.filecontent = self.filecontent:gsub('\n\n\n', '\n\n')
end

function File:remove_example()
  local prefix = self.basename:gsub("%.[ch]", "_")
  self:remove_function(prefix .. "example")
end

function File:remove_unit_test()
  self:remove_include("ensure.h")
  local prefix = self.basename:gsub("%.[ch]", "_")
  self:remove_function(prefix .. "unit_test")
end

function File:fa2xansi()
  self.filecontent = self.filecontent:gsub('#include "fa.h"',
                                           '#include "xansi.h"')
  self.filecontent = self.filecontent:gsub('fa_xfopen', 'xfopen')
  self.filecontent = self.filecontent:gsub('fa_xfclose', 'xfclose')
end

function File:ma2xansi()
  self.filecontent = self.filecontent:gsub('#include "ma.h"',
                                           '#include "xansi.h"')
  self.filecontent = self.filecontent:gsub('ma_malloc', 'xmalloc')
  self.filecontent = self.filecontent:gsub('ma_calloc', 'xcalloc')
  self.filecontent = self.filecontent:gsub('ma_realloc', 'xrealloc')
  self.filecontent = self.filecontent:gsub('cstr_dup', 'xstrdup')
  self.filecontent = self.filecontent:gsub('ma_free', 'free')
end

function File:simple_bioseq()
  self.filecontent = self.filecontent:gsub('#include "bioseq.h"',
                                           '#include "simple_bioseq.h"')
  self.filecontent = self.filecontent:gsub("Bioseq", "SimpleBioseq")
  self.filecontent = self.filecontent:gsub("bioseq_", "simple_bioseq_")
end

function File:add_test(test)
  assert(test)
  self.filecontent = self.filecontent ..
                     string.format("\n\n.PHONY: test\ntest: %s\n\t%s",
                                   self.progname, test)
end

function File:write(dir)
  assert(dir)
  local outfile = io.open(dir .. "/" .. self.basename, "w")
  outfile:write(self.filecontent)
  outfile:close()
end
