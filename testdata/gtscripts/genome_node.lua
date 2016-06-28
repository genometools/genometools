--[[
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2014      Sascha Steinbiss <sascha@steinbiss.name>
  Copyright (c) 2016      Daniel Standage <daniel.standage@gmail.com>
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

-- testing the Lua bindings for the GenomeNode interface

function count_children(parent)
  count = 0
  gfi = gt.feature_node_iterator_new(parent)
  curnode = gfi:next()
  while not(curnode == nil) do
      count = count + 1
      curnode = gfi:next()
  end
  return count
end

function table.contains(tab, element)
  for _, value in pairs(tab) do
    if value == element then
      return true
    end
  end
  return false
end

-- testing gt.feature_node_new
range = gt.range_new(1, 100)
rval, err = pcall(gt.feature_node_new, nil, nil, range:get_start(), range:get_end(), "+")
assert(not rval)
rval, err = pcall(gt.feature_node_new, "seqid", nil, range:get_start(), range:get_end(), "+")
assert(not rval)
rval, err = pcall(gt.feature_node_new, "seqid", "gene", "test", "+")
assert(not rval)
rval, err = pcall(gt.feature_node_new, "seqid", "gene", range:get_start(), range:get_end(), "plus")
assert(not rval)
assert(string.find(err, "strand string must have length 1"))
rval, err = pcall(gt.feature_node_new, "seqid", "gene", range:get_start(), range:get_end(), "p")
assert(not rval)
assert(string.find(err, "invalid strand"))
gn = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
assert(not gn:is_marked())
gn:mark()
assert(gn:is_marked())

parent = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
child  = gt.feature_node_new("seqid", "exon", range:get_start(), range:get_end(), "+")
parent:add_child(child)
assert(not parent:is_marked(parent))
assert(not parent:contains_marked(parent))
child:mark()
child  = nil; collectgarbage() -- being nasty
assert(not parent:is_marked(parent))
assert(parent:contains_marked(parent))

-- testing genome_node:get_filename
rval, fn = pcall(gn.get_filename, gn)
assert(rval)
assert(string.find(fn, "^generated$"))

-- testing feature_node:add_child
parent = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
child  = gt.feature_node_new("seqid2", "exon", range:get_start(), range:get_end(), "+")
rval, err = pcall(GenomeTools_genome_node.add_child, parent, child)
assert(not rval)
assert(string.find(err, "sequence ID different"))

-- testing genome_node:remove_leaf
-- testing removal of leaves which are direct children
parent = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
child  = gt.feature_node_new("seqid", "exon", range:get_start(), range:get_end(), "+")
parent:add_child(child)
child  = gt.feature_node_new("seqid", "exon", range:get_start(), range:get_end(), "+")
parent:add_child(child)
assert(count_children(parent) == 3)
parent:remove_leaf(child)
assert(count_children(parent) == 2)
parent:add_child(child)
assert(count_children(parent) == 3)
-- testing removal of leaves which are non-direct children
newchild = gt.feature_node_new("seqid", "exon", range:get_start(), range:get_end(), "+")
child:add_child(newchild)
assert(count_children(parent) == 4)
parent:remove_leaf(newchild)
assert(count_children(parent) == 3)

-- testing get_children
parent = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
child  = gt.feature_node_new("seqid", "exon", range:get_start(), range:get_end(), "+")
parent:add_child(child)
child2  = gt.feature_node_new("seqid", "exon", range:get_start()+1, range:get_end(), "+")
parent:add_child(child2)
out = {}
for i in parent:get_children() do
  table.insert(out, i)
end
assert(#out == 3)
assert(out[1] == parent)
assert(out[2] == child)
assert(out[3] == child2)
out = {}
for i in child2:get_children() do
  table.insert(out, i)
end
assert(#out == 1)
assert(out[3] == child1)
node = gt.sequence_node_new("foo", "acgatcgatcgcat")
rval, err = pcall(GenomeTools_genome_node.get_children, node)
assert(not rval)
assert(string.find(err, "not a feature node"))
node = gt.comment_node_new("foobar")
rval, err = pcall(GenomeTools_genome_node.get_children, node)
assert(not rval)
assert(string.find(err, "not a feature node"))
node = gt.region_node_new("foo", 1, 200)
rval, err = pcall(GenomeTools_genome_node.get_children, node)
assert(not rval)
assert(string.find(err, "not a feature node"))

-- testing get_direct_children()
child3  = gt.feature_node_new("seqid", "foo", range:get_start()+1, range:get_end(), "+")
child:add_child(child3)
out = {}
for i in parent:get_direct_children() do
  table.insert(out, i)
end
assert(#out == 2)
assert(out[1] == child)
assert(out[2] == child2)
out = {}
for i in child:get_direct_children() do
  table.insert(out, i)
end
assert(#out == 1)
assert(out[1] == child3)
out = {}
for i in child2:get_direct_children() do
  table.insert(out, i)
end
assert(#out == 0)
node = gt.sequence_node_new("foo", "acgatcgatcgcat")
rval, err = pcall(GenomeTools_genome_node.get_direct_children, node)
assert(not rval)
assert(string.find(err, "not a feature node"))
node = gt.comment_node_new("foobar")
rval, err = pcall(GenomeTools_genome_node.get_direct_children, node)
assert(not rval)
assert(string.find(err, "not a feature node"))
node = gt.region_node_new("foo", 1, 200)
rval, err = pcall(GenomeTools_genome_node.get_direct_children, node)
assert(not rval)
assert(string.find(err, "not a feature node"))

-- testing set_range()
rng2 = gt.range_new(2, 200)
node = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
assert(node:get_range():get_start() == 1)
assert(node:get_range():get_end() == 100)
node:set_range(rng2)
assert(node:get_range():get_start() == 2)
assert(node:get_range():get_end() == 200)
rval, err = pcall(GenomeTools_genome_node.set_range, node, "foo")
assert(not rval)
assert(string.find(err, "bad argument"))
rval, err = pcall(GenomeTools_genome_node.set_range, node, nil)
assert(not rval)
assert(string.find(err, "bad argument"))
node:set_range(range:join(rng2))
assert(node:get_range():get_start() == 1)
assert(node:get_range():get_end() == 200)
node = gt.sequence_node_new("foo", "acgatcgatcgcat")
node:set_range(rng2)
assert(node:get_range():get_start() == 0)
assert(node:get_range():get_end() == 0)
node = gt.comment_node_new("foobar")
node:set_range(rng2)
assert(node:get_range():get_start() == 0)
assert(node:get_range():get_end() == 0)
node = gt.region_node_new("foo", 1, 200)
node:set_range(rng2)
assert(node:get_range():get_start() == 2)
assert(node:get_range():get_end() == 200)

-- testing (get,add,set,remove)_attribute and attribute_pairs()

node = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
out = {}
n = 0
for k,v in node:attribute_pairs() do
  out[k] = v
  n = n + 1
end
assert(n == 0)
assert(not node:get_attribute("test"))
node:add_attribute("test","foo")
assert(node:get_attribute("test") == "foo")
out = {}
n = 0
for k,v in node:attribute_pairs() do
  out[k] = v
  n = n + 1
end
assert(n == 1)
assert(out.test == "foo")
rval, err = pcall(GenomeTools_genome_node.add_attribute, node, "test", "foo")
assert(not rval)
assert(string.find(err, "already present"))
node:set_attribute("test","bar")
assert(node:get_attribute("test") == "bar")
node:set_attribute("test", "baz")
assert(node:get_attribute("test") == "baz")
node:set_attribute("bar", "baz")
assert(node:get_attribute("bar") == "baz")
out = {}
n = 0
for k,v in node:attribute_pairs() do
  out[k] = v
  n = n + 1
end
assert(n == 2)
assert(out.test == "baz")
assert(out.bar == "baz")
rval, err = pcall(GenomeTools_genome_node.remove_attribute, node, "qqq")
assert(not rval)
assert(string.find(err, "not present"))
node:remove_attribute("test")
assert(node:get_attribute("test") == nil)
out = {}
n = 0
for k,v in node:attribute_pairs() do
  out[k] = v
  n = n + 1
end
assert(n == 1)
assert(out.test == nil)
assert(out.bar == "baz")
node:remove_attribute("bar")
assert(node:get_attribute("bar") == nil)
out = {}
n = 0
for k,v in node:attribute_pairs() do
  out[k] = v
  n = n + 1
end
assert(n == 0)
assert(out.test == nil)
assert(out.bar == nil)
node:add_attribute("test","foo")
assert(node:get_attribute("test") == "foo")
out = {}
n = 0
for k,v in node:attribute_pairs() do
  out[k] = v
  n = n + 1
end
assert(n == 1)
assert(out.test == "foo")

-- testing has_child_of_type
assert(parent:has_child_of_type("exon"))
assert(not parent:has_child_of_type("gene"))
assert(not parent:has_child_of_type("intron"))
assert(not child:has_child_of_type("gene"))
assert(not child:has_child_of_type("exon"))

-- testing count_children_of_type
range = gt.range_new(1000, 2000)
parent = gt.feature_node_new("seqid", "mRNA", 1000, 2000, "-")
child  = gt.feature_node_new("seqid", "exon", 1000, 1200, "-")
parent:add_child(child)
child  = gt.feature_node_new("seqid", "exon", 1350, 1800, "-")
parent:add_child(child)
child  = gt.feature_node_new("seqid", "exon", 1980, 1200, "-")
parent:add_child(child)
child  = gt.feature_node_new("seqid", "intron", 1201, 1349, "-")
parent:add_child(child)
child  = gt.feature_node_new("seqid", "intron", 1801, 1979, "-")
parent:add_child(child)
assert(parent:count_children_of_type("exon") == 3)
assert(parent:count_children_of_type("intron") == 2)
assert(parent:count_children_of_type("b0gu$") == 0)

-- testing change_seqid
node = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
assert(node:get_seqid() == "seqid")
node:change_seqid("foo")
assert(node:get_seqid() == "foo")
node:change_seqid("seqid")
assert(node:get_seqid() == "seqid")

-- testing set_type
node = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
assert(node:get_type() == "gene")
node:set_type("pseudogene")
assert(node:get_type() == "pseudogene")
node:set_type(1)
assert(node:get_type() == "1")
node:set_type("gene")
assert(node:get_type() == "gene")
rval, err = pcall(GenomeTools_genome_node.set_type, node, nil)
assert(not rval)
assert(err:match("string expected, got nil"))
assert(node:get_type() == "gene")
rval, err = pcall(GenomeTools_genome_node.set_type, node, "")
assert(not rval)
assert(err:match("must not be empty"))
assert(node:get_type() == "gene")

-- testing set_strand
node = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
assert(node:get_strand() == "+")
node:set_strand("-")
assert(node:get_strand() == "-")
node:set_strand(".")
assert(node:get_strand() == ".")
node:set_strand("?")
assert(node:get_strand() == "?")
node:set_strand("+")
assert(node:get_strand() == "+")
node:set_strand("-")
assert(node:get_strand() == "-")
rval, err = pcall(GenomeTools_genome_node.set_strand, node, "X")
assert(not rval)
assert(err:match("must be one of"))
assert(node:get_strand() == "-")
rval, err = pcall(GenomeTools_genome_node.set_strand, node, "+-")
assert(not rval)
assert(err:match("must be one of"))
assert(node:get_strand() == "-")

-- testing set_score
node = gt.feature_node_new("seqid", "gene", range:get_start(), range:get_end(), "+")
assert(node:get_score() == nil)
node:set_score(1)
assert(node:get_score() == 1)
rval, err = pcall(GenomeTools_genome_node.set_score, node, "foo")
assert(not rval)
assert(node:get_score() == 1)
rval, err = pcall(GenomeTools_genome_node.set_score, node, nil)
assert(not rval)
assert(node:get_score() == 1)

-- testing set_phase
node = gt.feature_node_new("seqid", "CDS", range:get_start(), range:get_end(), "+")
assert(node:get_phase() == ".")
node:set_phase(1)
assert(node:get_phase() == "1")
node:set_phase(0)
assert(node:get_phase() == "0")
node:set_phase(2)
assert(node:get_phase() == "2")
node:set_phase(".")
assert(node:get_phase() == ".")
node:set_phase(1)
assert(node:get_phase() == "1")
rval, err = pcall(GenomeTools_genome_node.set_phase, node, "foo")
assert(not rval)
assert(node:get_phase() == "1")
rval, err = pcall(GenomeTools_genome_node.set_phase, node, nil)
assert(not rval)
assert(node:get_phase() == "1")

-- testing gt.region_node_new
range = gt.range_new(1, 100)
rval, err = pcall(gt.region_node_new, nil, range:get_start(), range:get_end())
assert(not rval)
rval, err = pcall(gt.region_node_new, "chr1", "test")
assert(not rval)
gn = gt.region_node_new("chr1", range:get_start(), range:get_end())

-- testing gt.meta_node_new
rval, err = pcall(gt.meta_node_new, nil, "test")
assert(not rval)
rval, err = pcall(gt.region_node_new, "foo", nil)
assert(not rval)
gn = gt.meta_node_new("foo","bar")
assert(gn:get_directive() == "foo")
assert(gn:get_data() == "bar")
gn = gt.meta_node_new("foo", nil)
assert(gn:get_directive() == "foo")
assert(not gn:get_data())

-- testing gt.comment_node_new
rval, err = pcall(gt.comment_node_new, nil)
assert(not rval)
cn = gt.comment_node_new("bar")
assert(cn:get_comment() == "bar")
cn = gt.comment_node_new(42)
assert(cn:get_comment() == "42")
cn:get_range()

-- testing gt.sequence_node_new
rval, err = pcall(gt.sequence_node_new, nil)
assert(not rval)
rval, err = pcall(gt.sequence_node_new, nil, "foo")
assert(not rval)
rval, err = pcall(gt.sequence_node_new, "foo", nil)
assert(not rval)
sn = gt.sequence_node_new("bar", "CTGA")
assert(sn:get_sequence() == "CTGA")
assert(sn:get_sequence_length() == 4)
assert(sn:get_description() == "bar")
sn = gt.sequence_node_new("bar", "")
assert(sn:get_sequence() == "")
assert(sn:get_sequence_length() == 0)
assert(sn:get_description() == "bar")
