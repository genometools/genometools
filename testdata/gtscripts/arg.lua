--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
]]

-- testing the arg table

io.write(string.format("arg[0]=%s\n", arg[0]))
for i,v in ipairs(arg) do
  io.write(string.format("arg[%d]=%s\n", i, v))
end
