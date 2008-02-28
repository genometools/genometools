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
require 'lfs'

function usage()
  io.stderr:write(string.format("Usage: %s PNG_dir reality_file " ..
                                "prediction_file\n", arg[0]))
  io.stderr:write("Evaluate prediction_file against reality_file and write " ..
                  "out PNGs to PNG_dir.\n")
  os.exit(1)
end

if #arg == 3 then
  png_dir  = arg[1]
  reality_file = arg[2]
  prediction_file = arg[3]
  -- make sure png_dir is a directory or create it
  rval, err = lfs.attributes(png_dir, "mode")
  if rval then
    if rval ~= "directory" then
      io.stderr:write(string.format("PNG_dir '%s' is not a directory\n",
                                    png_dir))
      os.exit(1)
    end
  else
    -- not successfull, try to create directory
    rval, err = lfs.mkdir(png_dir)
    if not rval then
      io.stderr:write(string.format("could not create directory '%s': %s",
                                    png_dir, err))
      os.exit(1)
    end
  end
else
  usage()
end

function write_marked_regions(seqid, filenumber, maxdist)
  assert(seqid)
  local marked = feature_index:get_marked_regions(seqid, maxdist)
  for _,range in ipairs(marked) do
    local filename = png_dir .. "/" .. filenumber .. ".png"
    io.write(string.format("writing file '%s'\n", filename))
    feature_index:render_to_png(seqid, range, filename, width)
    filenumber = filenumber + 1
  end
  return filenumber
end

-- process input files
reality_stream = gt.gff3_in_stream_new_sorted(reality_file)
prediction_stream = gt.gff3_in_stream_new_sorted(prediction_file)
stream_evaluator = gt.stream_evaluator_new(reality_stream, prediction_stream)

feature_index = gt.feature_index_new()
feature_visitor = gt.feature_visitor_new(feature_index)

stream_evaluator:evaluate(feature_visitor)
stream_evaluator:show()

-- write results
filenumber = 1
width = 1600
for _, seqid in ipairs(feature_index:get_seqids()) do
  print(string.format("seqid '%s'", seqid))
  filenumber = write_marked_regions(seqid, filenumber)
end

-- get ready for interactive mode
fi = feature_index
features = fi:get_all_features()
marked_features = gt.features_get_marked(features)
gt.export()
