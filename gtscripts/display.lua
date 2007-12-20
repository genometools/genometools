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

require 'gtlua'

function usage()
  io.stderr:write(string.format("Usage: %s PNG_dir\n", arg[0]))
  io.stderr:write("Display numbered PNG files in PNG_dir.\n")
  os.exit(1)
end

if #arg == 1 then
  png_dir  = arg[1]
  -- make sure png_dir is a directory
  rval, err = lfs.attributes(png_dir, "mode")
  if rval ~= "directory" then
    io.stderr:write(string.format("PNG_dir '%s' is not a directory\n", png_dir))
    os.exit(1)
  end
else
  usage()
end

filenumber = 1
while (true) do
  local filename = png_dir .. "/" .. filenumber .. ".png"
  if gt.file_exists(filename) then
    gt.display(filename)
    filenumber = filenumber + 1
  else
    break
  end
end
