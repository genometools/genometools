#
# Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2009 Center for Bioinformatics, University of Hamburg
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

# test CommentNodes
cn  = GT::CommentNode.create("Comment")
cn2 = GT::CommentNode.create(333)

unless cn.get_comment == "Comment"
  raise(TestFailedError, "CommentNode.get_comment")
end
unless cn2.get_comment == "333"
  raise(TestFailedError, "CommentNode.get_comment")
end

# test SequenceNodes
SEQUENCE = 'AGATATAGA'
DESC = 'testdesc'
begin
  sn = GT::SequenceNode.create("foo")
rescue ArgumentError
  # expected
else
  raise(TestFailedError, "no failure on not enough arguments")
end
sn = GT::SequenceNode.create(DESC, SEQUENCE)
unless sn.get_description == 'testdesc'
  raise(TestFailedError, "SequenceNode.get_description")
end
unless (seq = sn.get_sequence) == 'AGATATAGA'
  raise(TestFailedError, "SequenceNode.get_sequence")
end
unless sn.get_sequence_length == 9
  raise(TestFailedError, "SequenceNode.get_sequence_length")
end
unless sn.get_sequence_length == seq.length
  raise(TestFailedError, "SequenceNode.get_sequence_length")
end

# test RegionNodes
begin
  rn = GT::RegionNode.create("foo", 100, 50)
rescue ArgumentError
  # expected
else
  raise(TestFailedError, "no failure on start > stop")
end
# does not have any exposed methods yet (see _api.h)

# test MetaNodes
cn  = GT::MetaNode.create("directive", "data")
cn2 = GT::MetaNode.create(333, 444)

unless cn.get_directive == "directive"
  raise(TestFailedError, "MetaNode.get_directive")
end
unless cn2.get_directive == "333"
  raise(TestFailedError, "CommentNode.get_directive")
end
unless cn.get_data == "data"
  raise(TestFailedError, "MetaNode.get_data")
end
unless cn2.get_data == "444"
  raise(TestFailedError, "CommentNode.get_data")
end

# test EOFNodes
rn = GT::EOFNode.create()
if rn.nil? then
  raise(TestFailedError, "EOFNode.create")
end
# does not have any exposed methods yet (see _api.h)
