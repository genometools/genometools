#!/usr/bin/env ruby
#
# Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg
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
  $:<< "."            # favor the local stest version
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
  $cur=File.join(Dir.pwd, "..", "")
end

$transdir=File.join(Dir.pwd, "..", "gtdata" , "trans", "")

$scriptsdir=File.join(Dir.pwd, "..", "scripts", "")

if $arguments["gtruby"] then
  $gtruby=File.join($arguments["gtruby"], "")
else
  $gtruby=File.join(Dir.pwd, "..", "gtruby", "")
end

if $arguments["gtpython"] then
  $gtpython=File.join($arguments["gtpython"], "")
else
  $gtpython=File.join(Dir.pwd, "..", "gtpython", "")
end

if $arguments["gttestdata"] then
  $gttestdata=File.join($arguments["gttestdata"], "")
end

$systemname=`uname -s`
$systemname.chomp!

$MEMCHECK_SUPPORTEDPLATFORMS = ["Linux", "Darwin"]

# define helper function
def run_test(str, opts = {})
  if $arguments["memcheck"] then
    if $MEMCHECK_SUPPORTEDPLATFORMS.include?($systemname) then
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

def run_python(str, opts = {})
  run("env PYTHONPATH=#{$gtpython} LD_LIBRARY_PATH=#{$cur}/lib python " + \
      "#{$path}#{str}", opts)
end

# include the actual test modules
require 'gt_bed_to_gff3_include'
require 'gt_cds_include'
require 'gt_chseqids_include'
require 'gt_consensus_sa_include'
require 'gt_csa_include'
require 'gt_csr_include.rb'
require 'gt_encseq_include'
require 'gt_eval_include'
require 'gt_extractfeat_include'
require 'gt_featureindex_include'
require 'gt_fingerprint_include'
require 'gt_genomediff_include'
require 'gt_gff3_include'
require 'gt_gff3validator_include'
require 'gt_gtf_to_gff3_include'
require 'gt_hop_include'
require 'gt_id_to_md5_include'
require 'gt_include'
require 'gt_interfeat_include'
require 'gt_ltrdigest_include'
require 'gt_ltrharvest_include'
require 'gt_magicmatch_include'
require 'gt_matchtool_include'
require 'gt_mathsupport_include'
require 'gt_md5_to_id_include'
require 'gt_merge_include'
require 'gt_mergefeat_include'
require 'gt_mgth_include'
require 'gt_mmapandread_include'
require 'gt_orffinder_include'
require 'gt_python_include'
require 'gt_readjoiner_include'
require 'gt_readreads_include'
require 'gt_regioncov_include'
require 'gt_ruby_include'
require 'gt_sambam_include'
require 'gt_script_filter_include'
require 'gt_scripts_include'
require 'gt_select_include'
require 'gt_seq_include'
require 'gt_seqbuffer_include'
require 'gt_seqfilter_include'
require 'gt_seqids_include'
require 'gt_seqlensort_include'
require 'gt_seqmutate_include'
require 'gt_seqorder_include'
require 'gt_seqstat_include'
require 'gt_seqtransform_include'
require 'gt_sequniq_include'
require 'gt_shredder_include'
require 'gt_simreads_include'
require 'gt_splicesiteinfo_include'
require 'gt_splitfasta_include'
require 'gt_stat_include'
require 'gt_uniq_include'
if not $arguments["nocairo"] then
  require 'gt_sketch_include'
end

#now the test cases for the tools implemented or supervised by SK.

require 'gt_chain2dim_include'
require 'gt_cmpprjfiles_include'
require 'gt_env_options_include'
require 'gt_extractseq_include'
require 'gt_idxsearch_include'
require 'gt_mergeesa_include'
require 'gt_packedindex_include'
require 'gt_sortbench_include'
require 'gt_suffixerator_include'
require 'gt_encseq2spm_include'
require 'gt_tallymer_include'
require 'gt_trieins_include'
require 'scripts_include'

if $arguments["gcov"] then
  require 'gcov_include' # must be last
end

#we now have all tests in $testsuite.

if $arguments["threads"] then
  $testsuite.nof_threads = $arguments["threads"].to_i
end
$testsuite.run
