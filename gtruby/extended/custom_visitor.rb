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

require 'dl/import'
require 'gthelper'
require 'extended/comment_node'
require 'extended/feature_node'
require 'extended/sequence_node'
require 'extended/region_node'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtNodeVisitor* gt_script_wrapper_visitor_new(void*, void*, void*,
                                                       void*, void*)"

  class CustomVisitor
    def visitor_func_generic(node_class, node_ptr, err_ptr, method_name)
      node = node_class.new(node_ptr, true)
      err = GT::Error.new(err_ptr)
      begin
        begin
          ret = self.send("visit_#{method_name}_node", node)
        rescue NoMethodError
          ret = 0
        end
      rescue Error => msg
        err.set(msg)
        ret = -1
      end
      ret.nil? ? 0 : ret.to_i
    end

    def initialize()
      @feature_node_cb = DL.callback("IPP") do |fn_ptr, err_ptr|
        self.visitor_func_generic(GT::FeatureNode, fn_ptr, err_ptr, "feature")
      end
      @comment_node_cb = DL.callback("IPP") do |cn_ptr, err_ptr|
        self.visitor_func_generic(GT::CommentNode, cn_ptr, err_ptr, "comment")
      end
      @region_node_cb = DL.callback("IPP") do |rn_ptr, err_ptr|
        self.visitor_func_generic(GT::RegionNode, rn_ptr, err_ptr, "region")
      end
      @sequence_node_cb = DL.callback("IPP") do |sn_ptr, err_ptr|
        self.visitor_func_generic(GT::CommentNode, sn_ptr, err_ptr, "sequence")
      end
      @genome_visitor = GT.gt_script_wrapper_visitor_new(@comment_node_cb,
                                                         @feature_node_cb,
                                                         @region_node_cb,
                                                         @sequence_node_cb,
                                                         nil)
      @genome_visitor.free = GT::symbol("gt_node_visitor_delete", "0P")
    end

    def to_ptr
      @genome_visitor
    end
  end
end
