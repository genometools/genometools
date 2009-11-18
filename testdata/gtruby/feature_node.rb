#
# Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008 Center for Bioinformatics, University of Hamburg
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


class TestFailedError < Exception
end

fn = GT::FeatureNode.create("test", "type", 100, 500, "+")
fn.add_attribute("test","testval")
fn.add_attribute("test2","testval2")
if fn.score_is_defined? then
  raise TestFailedError
end
fn.set_score(2)
if !fn.score_is_defined? then
  raise TestFailedError
end
if fn.get_score != 2 then
  raise TestFailedError
end
fn.unset_score()
if fn.score_is_defined? then
  raise TestFailedError
end
if fn.has_type?("foo") then
  raise TestFailedError
end
if !fn.has_type?("type") then
  raise TestFailedError
end
if !fn.get_strand == "+" then
  raise TestFailedError
end
fn2 = GT::FeatureNode.create("test", "type2", 200,300,"+")
fn.add_child(fn2)
num_attrs = 0
fn.each_attribute do |tag, val|
  if val != fn.get_attribute(tag) then
    raise TestFailedError
  end
  num_attrs += 1
end
if !num_attrs == 2 then
  raise TestFailedError
end
if fn.get_phase != 3 then   #undefined
  raise TestFailedError
end
fn.set_phase(0)
if fn.get_phase != 0 then   #zero
  raise TestFailedError
end
if fn.get_filename != "generated" then
  raise TestFailedError
end

begin
  fn.add_child(GT::FeatureNode.create("nottest", "foo", 100, 200, "+"))
rescue GT::GTError => msg
  # expect exception
else
  raise TestFailedError
end

fni = GT::FeatureNodeIteratorDepthFirst.new(fn)
num_features = 0
tfn = fni.next
while tfn do
  num_features += 1
  tfn = fni.next
end
if num_features != 2 then
  raise TestFailedError
end

fn3 = GT::FeatureNode.create("test", "type3", 250, 300,"+")
fn.add_child(fn3)
fni = GT::FeatureNodeIteratorDepthFirst.new(fn)
num_features = 0
tfn = fni.next
while tfn do
  num_features += 1
  tfn = fni.next
end
if num_features != 3 then
  raise TestFailedError
end

types = []
fn.traverse_dfs do |tfn|
  types.push(tfn.get_type())
end
if types != ["type", "type2", "type3"] then
  raise TestFailedError
end

types = []
fn.traverse_direct do |tfn|
  types.push(tfn.get_type())
end
if types != ["type2", "type3"] then
  raise TestFailedError
end
