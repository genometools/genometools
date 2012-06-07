--[[
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

require 'gtlua'

function usage()
  io.stderr:write(string.format("Usage: %s <positive integer>\n" , arg[0]))
  os.exit(1)
end

function outputmulticseq(indent, linewidth, times)
  for i = 1, times-1 do
    io.write('c')
    if indent >= linewidth then
      io.write('\n')
      indent = 0
    else
      indent = indent + 1
    end
  end
  return indent
end

function outputrepetitivesequence(iter)
  indent = 1;
  linewidth = 70;
  midseq = "acaccaccc"

  io.write(string.format(">%s for i=%d\na","ac^{i^2}acac^2ac^3 ... ac^ia",iter))
  indent = outputmulticseq(indent,linewidth,iter * iter)
  assert (indent < linewidth)
  if indent + string.len(midseq) > linewidth then
    rest = linewidth - indent
    io.write(string.format("%*.*s\n",rest,rest,midseq))
    if rest < string.len(midseq) then
      io.write(string.format("%s",midseq + rest))
    end
    indent = 0
  else
    io.write(string.format("acaccaccc"))
    indent = indent + string.len(midseq)
  end
  indent = outputmulticseq(indent,linewidth,iter)
  io.write("\n")
end

if #arg == 1 then
  i = arg[1]
  outputrepetitivesequence(i)
else
  usage()
end

gt.export()
