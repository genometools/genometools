#
# Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2010 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

require 'gtruby'

def assert_raises
  begin
    yield
    raise
  rescue GT::GTError, ArgumentError => e
    #pass
  end
end

assert_raises do
  GT::Range.new()
end

assert_raises do
  GT::Range.new(1)
end

assert_raises do
  GT::Range.new(1, nil)
end

assert_raises do
  GT::Range.new(nil, 1)
end

assert_raises do
  GT::Range.new(20, 1)
end

assert_raises do
  GT::Range.new(-1, 1)
end

assert_raises do
  GT::Range.new(1, -1)
end

assert_raises do
  GT::Range.new(-1, -1)
end

assert_raises do
  GT::Range.new(-1, -2)
end

raise unless !GT::Range.new(1, 2).nil?
raise unless !GT::Range.new(1, 1).nil?

rng = GT::Range.new(100, 200)

assert_raises do
  rng.start = 300
end

assert_raises do
  rng.start = -100
end

assert_raises do
  rng.end = 50
end

assert_raises do
  rng.end = -350
end

rng.start = 50
raise unless rng.start == 50
rng.start = 150
raise unless rng.start == 150

rng.end = 250
raise unless rng.end == 250
rng.end = 151
raise unless rng.end == 151
rng.end = 150
raise unless rng.end == 150

assert_raises do
  rng.set(1)
end

assert_raises do
  rng.set(1, nil)
end

assert_raises do
  rng.set(nil, 1)
end

assert_raises do
  rng.set(20, 1)
end

assert_raises do
  rng.set(-1, 1)
end

assert_raises do
  rng.set(1, -1)
end

assert_raises do
  rng.set(-1, -1)
end

assert_raises do
  rng.set(-1, -2)
end
