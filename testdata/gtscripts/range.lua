--[[
  Copyright (c) 2007 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  Copyright (c) 2014 Genome Research Ltd.

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

-- testing the Lua bindings for the Range class

range = gt.range_new(1, 1000)
assert(range:get_start() == 1)
assert(range:get_end() == 1000)

rval, err = pcall(gt.range_new, 1000, 1)
assert(not rval)
assert(string.find(err, "must be <= endpos"))

ranges = {}
for i = 1, 100 do
  range = gt.range_new(i, i+1)
  table.insert(ranges, range)
end
ranges = gt.ranges_sort(ranges)
assert(gt.ranges_are_sorted(ranges))

-- join

range_a = gt.range_new(1, 1000)
range_b = gt.range_new(400, 3000)
range_c = range_a:join(range_b)
assert(range_c:get_start() == 1)
assert(range_c:get_end() == 3000)

range_b = gt.range_new(1, 1000)
range_a = gt.range_new(2000, 3000)
range_c = range_a:join(range_b)
assert(range_c:get_start() == 1)
assert(range_c:get_end() == 3000)

-- contains

range_a = gt.range_new(1, 1000)
range_b = gt.range_new(1, 300)
assert(range_a:contains(range_b))
assert(not range_b:contains(range_a))

-- within

range = gt.range_new(1, 1000)
assert(range_a:within(300))
assert(not range_a:within(1300))

-- overlaps

range_a = gt.range_new(1, 1000)
range_b = gt.range_new(400, 3000)
range_c = gt.range_new(2000, 3000)
assert(range_a:overlap(range_b))
assert(range_b:overlap(range_a))
assert(not range_a:overlap(range_c))

-- length

range = gt.range_new(1, 1000)
assert(range:length() == 1000)
range = gt.range_new(1, 1)
assert(range:length() == 1)