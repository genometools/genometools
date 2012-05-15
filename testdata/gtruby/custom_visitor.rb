#
# Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2012 Center for Bioinformatics, University of Hamburg
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

module GT
  class TestVisitor < CustomVisitor
    attr_reader :sn, :en, :rn, :cn, :fn, :mn

    def initialize
      super()
      @sn = nil
      @en = nil
      @rn = nil
      @cn = nil
      @fn = nil
      @mn = nil
    end

    def visit_feature_node(fn)
      @fn = fn
      0
    end

    def visit_sequence_node(sn)
      @sn = sn
      0
    end

    def visit_comment_node(cn)
      @cn = cn
      0
    end

    def visit_region_node(rn)
      @rn = rn
      0
    end

    def visit_meta_node(mn)
      @mn = mn
      0
    end

    def visit_eof_node(en)
      @en = en
      0
    end
  end
end

fn = GT::FeatureNode.create("foo", "gene", 100, 10000, "+")
cn = GT::CommentNode.create("comment")
rn = GT::RegionNode.create("foo", 100, 2000)
sn = GT::SequenceNode.create("foo", "AGATATAGA")
en = GT::EOFNode.create
mn = GT::MetaNode.create("foo", "bar")

v = GT::TestVisitor.new

raise unless v.sn.nil?
sn.accept(v)
raise if v.sn.nil? or v.sn != sn

raise unless v.cn.nil?
cn.accept(v)
raise if v.cn.nil? or v.cn != cn

raise unless v.fn.nil?
fn.accept(v)
raise if v.fn.nil? or v.fn != fn

raise unless v.en.nil?
en.accept(v)
raise if v.en.nil? or v.en != en

raise unless v.rn.nil?
rn.accept(v)
raise if v.rn.nil? or v.rn != rn

raise unless v.mn.nil?
mn.accept(v)
raise if v.mn.nil? or v.mn != mn
