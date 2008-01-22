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

-- Returns true if the given array of <features> contains a marked feature,
-- false otherwise.
function features_contain_marked(features)
  assert(features)
  for _, feature in ipairs(features) do
    if feature:contains_marked() then
      return true
    end
  end
  return false
end

-- Print the given array of <features> to stdout.
function features_show(features)
  assert(features)
  local gff3_visitor = gt.gff3_visitor_new()
  for _, features in ipairs(features) do
    features:show(gff3_visitor)
  end
end

-- Print all marked <features> (an array) to stdout.
function features_show_marked(features)
  assert(features)
  if features_contain_marked(features) then
    for _, feature in ipairs(features) do
      feature:show_marked()
    end
  end
end
