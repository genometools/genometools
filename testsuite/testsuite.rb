#!/usr/bin/env ruby
#
# Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
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

# the GenomeTools test suite (employs ``stest'').

if $0 == __FILE__
  require 'stest'
  at_exit do
    OnError do exit 1 end
  end
end

# set some global variables
if $arguments["path"] then
  $path=File.join($arguments["path"], "")
else
  $path=""
end

if $arguments["testdata"] then
  $testdata=File.join($arguments["testdata"], "")
else
  $testdata=File.join(Dir.pwd, "..", "testdata", "")
end

if $arguments["bin"] then
  $bin=File.join($arguments["bin"], "")
else
  $bin=File.join(Dir.pwd, "..", "bin", "")
end

if $arguments["cur"] then
  $cur=$arguments["cur"]
else
  $cur=File.join(Dir.pwd, "..")
end

$transdir=File.join(Dir.pwd, "..", "gtdata" , "trans", "")

$scriptsdir=File.join(Dir.pwd, "..", "scripts", "")

if $arguments["gtruby"] then
  $gtruby=File.join($arguments["gtruby"], "")
else
  $gtruby=File.join(Dir.pwd, "..", "gtruby", "")
end

if $arguments["gttestdata"] then
  $gttestdata=File.join($arguments["gttestdata"], "")
end

$systemname=`uname -s`
$systemname.chomp!

# define helper function
def run_test(str, opts = {})
  if $arguments["memcheck"] then
    if $systemname == "Linux" then
      $memcheck = "valgrind --tool=memcheck --suppressions="+
                  File.join($testdata, "gt.supp")+
                  " --leak-check=yes --error-exitcode=1 -q"
    elsif $systemname == "OpenBSD" then
      $memcheck = "env MALLOC_OPTIONS='GJ'"
    end
  else
    $memcheck = ""
  end
  run("#{$memcheck} #{$path}#{str}", opts)
end

def run_ruby(str, opts = {})
  run("env LD_LIBRARY_PATH=#{$cur}/lib ruby -I #{$gtruby} #{$path}#{str}", opts)
end

# include the actual test modules
require 'gt_include'
require 'gt_bioseq_include'
require 'gt_chseqids_include'
require 'gt_cds_include'
require 'gt_csa_include'
require 'gt_eval_include'
require 'gt_exercise_include'
require 'gt_extractfeat_include'
require 'gt_extractseq_include'
require 'gt_filter_include'
require 'gt_fingerprint_include'
require 'gt_gff3_include'
require 'gt_gtf_to_gff3_include'
require 'gt_ltrharvest_include'
require 'gt_magicmatch_include'
require 'gt_merge_include'
require 'gt_mgth_include'
require 'gt_mmapandread_include'
require 'gt_mutate_include'
require 'gt_packedindex_include'
require 'gt_regioncov_include'
require 'gt_ruby_include'
require 'gt_scripts_include'
require 'gt_seqiterator_include'
require 'gt_sequniq_include'
require 'gt_shredder_include'
require 'gt_splicesiteinfo_include'
require 'gt_splitfasta_include'
require 'gt_stat_include'
require 'gt_suffixerator_include'
require 'gt_idxsearch_include'
require 'gt_trieins_include'
require 'gt_uniq_include'
if $arguments["libgtview"] then
  require 'gt_view_include'
end
require 'gt_env_options_include'
require 'scripts_include'

if $arguments["gcov"] then
  require 'gcov_include' # must be last
end
