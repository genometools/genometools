#!/usr/bin/env ruby
# 
# Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2006 Center for Bioinformatics, University of Hamburg 
# See LICENSE file or http://genometools.org/license.html for license details.
# 

# the GenomeThreader test suite (employs ``stest'').

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

$systemname=`uname -s`
$systemname.chomp!

# define helper function
def run_test(str, opts = {})
  if $arguments["memcheck"] then
    if $systemname == "Linux" then 
      memcheck = "valgrind --tool=memcheck --suppressions="+
                 File.join($testdata, "gt.supp")+
                 " --leak-check=yes -q"
    elsif $systemname == "OpenBSD" then
      memcheck = "env MALLOC_OPTIONS='GJ'"
    end
  else
    memcheck = ""
  end
  run("#{memcheck} #{$path}#{str}", opts)
end

# include the actual test modules
require 'gt_bioseq_include'
require 'gt_cds_include'
require 'gt_csa_include'
require 'gt_eval_include'
require 'gt_exercise_include'
require 'gt_extractfeat_include'
require 'gt_gff3_include'
require 'gt_gtf2gff3_include'
require 'gt_merge_include'
require 'gt_mmapandread_include'
