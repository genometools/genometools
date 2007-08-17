--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
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

-- testing bittab:and
