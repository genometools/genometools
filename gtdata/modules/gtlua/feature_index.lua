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

function GenomeTools_feature_index_mt:get_coverage(seqid, maxdist)
  assert(seqid)
  local maxdist = maxdist or 0
  local features = self:get_features_for_seqid(seqid)
  local starpos, endpos
  local minstartpos = nil
  local maxendpos = nil
  local ranges = {}
  local coverage = {}

  -- collect all feature ranges
  for i, feature in ipairs(features) do
    table.insert(ranges, feature:get_range())
  end
  -- sort feature ranges
  ranges = gt.ranges_sort(ranges)

  -- compute and store coverage
  for i, range in ipairs(ranges) do
    startpos, endpos = range:get_start(), range:get_end()
    if i == 1 then
      minstartpos = startpos
      maxendpos   = endpos
    else
      -- assert(startpos >= minstartpos)
      if (startpos > maxendpos + maxdist) then
        -- new region started
        table.insert(coverage, gt.range_new(minstartpos, maxendpos))
        minstartpos = startpos
        maxendpos   = endpos
      else
        -- continue old region
        maxendpos = (endpos > maxendpos) and endpos or maxendpos
      end
    end
  end
  -- add last region
  coverage[#coverage+1] = gt.range_new(minstartpos, maxendpos)
  return coverage
end

-- render PNG file <png_file> for <seqid> in <range> with optional <width>
function GenomeTools_feature_index_mt:render_to_png(png_file, seqid, range, width)
  assert(self and png_file and seqid and range)
  if not width then width = 1600 end
  local diagram = gt.diagram_new(self, seqid, range)
  local render =  gt.render_new()
  render:to_png(diagram, png_file, width)
end

-- show all sequence ids
function GenomeTools_feature_index_mt:show_seqids()
  for _,seqid in ipairs(feature_index:get_seqids()) do
    print(seqid)
  end
end
