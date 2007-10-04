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

-- testing the Lua bindings for the Bittab class

-- testing gt.bittab_new
b = gt.bittab_new(10)
rval, err = pcall(gt.bittab_new, 0)
assert(not rval)
assert(string.find(err, "must be > 0"))

-- testing bittab:set_bit
b:set_bit(1)
rval, err = pcall(b.set_bit, b, 10)
assert(not rval)
assert(string.find(err, "bit number too large"))
rval, err = pcall(b.set_bit, a, 1)
assert(not rval)
assert(string.find(err, "bittab expected"))

-- testing bittab:unset_bit
b:unset_bit(1);
rval, err = pcall(b.unset_bit, b, 10)
assert(not rval)
assert(string.find(err, "bit number too large"))

-- testing bittab:complement
src = gt.bittab_new(10)
src:set_bit(5)
src:set_bit(7)
prob = gt.bittab_new(11)
b:complement(src)
rval, err = pcall(b.complement, b, prob)
assert(not rval)
assert(string.find(err, "bittabs have different sizes"))

-- testing bittab:equal
b:equal(src)
rval, err = pcall(b.equal, b, prob)
assert(not rval)
assert(string.find(err, "bittabs have different sizes"))

-- testing bittab:and_equal and bittab:bit_is_set
a = gt.bittab_new(100)
b = gt.bittab_new(100)
a:set_bit(0)
a:set_bit(50)
b:set_bit(50)
b:set_bit(99)
a:and_equal(b)
assert(not a:bit_is_set(0))
assert(a:bit_is_set(50))
assert(not a:bit_is_set(99))
