#!/usr/bin/env ruby
#
# Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
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

if $arguments["seed"] then
  $SEED = $arguments["seed"]
else
  $SEED = rand(2**31)
end

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
  run("env G_DEBUG=gc-friendly G_SLICE=always-malloc GT_SEED=#{$SEED} " + \
      "#{$memcheck} #{$path}#{str}", opts)
end

def run_ruby(str, opts = {})
  run("env G_DEBUG=gc-friendly G_SLICE=always-malloc " + \
      "DYLD_LIBRARY_PATH=#{$cur}/lib LD_LIBRARY_PATH=#{$cur}/lib " + \
      "ruby -I #{$gtruby} #{$path}#{str}", opts)
end

def run_python(str, opts = {})
  run("env G_DEBUG=gc-friendly G_SLICE=always-malloc " + \
      "PYTHONPATH=#{$gtpython} DYLD_LIBRARY_PATH=#{$cur}/lib " + \
      "LD_LIBRARY_PATH=#{$cur}/lib python " + \
      "#{$path}#{str}", opts)
end

def with_environment(variables={})
  if block_given?
    old_values = variables.map{ |k,v| [k,ENV[k]] }
    begin
       variables.each{ |k,v| ENV[k] = v }
       result = yield
    ensure
      old_values.each{ |k,v| ENV[k] = v }
    end
    result
  else
    variables.each{ |k,v| ENV[k] = v }
  end
end

def python_tests_runnable?
  if not ENV.has_key?("PATH")
    return false
  end
  return false unless
    ENV['PATH'].split(File::PATH_SEPARATOR).any? do |directory|
      File.executable?(File.join(directory, "python"))
    end
  require "open3"
  runline = "python #{$gtpython}/gt/dlload.py"
  with_environment({"PYTHONPATH" => $gtpython, \
                    "DYLD_LIBRARY_PATH" => "#{$cur}/lib", \
                    "LD_LIBRARY_PATH" => "#{$cur}/lib"}) do
    stdin, stdout, stderr = Open3.popen3(runline)
    result = stderr.read.match(/OSError/).nil?
    result
  end
end

def ruby_tests_runnable?
  rv = RUBY_VERSION.match(/^(\d+\.\d+)/)
  if rv.nil? or rv[1].to_f >= 1.9 then
    return false
  end
  require "open3"
  runline = "ruby -I #{$gtruby} #{$gtruby}/gtdlload.rb"
  with_environment({"GTRUBY" => $gtruby, \
                    "DYLD_LIBRARY_PATH" => "#{$cur}/lib", \
                    "LD_LIBRARY_PATH" => "#{$cur}/lib"}) do
    stdin, stdout, stderr = Open3.popen3(runline)
    result = stderr.read.match(/RuntimeError/).nil?
    result
  end
end

# include the actual test modules
require 'gt_bed_to_gff3_include'
require 'gt_cds_include'
require 'gt_chseqids_include'
require 'gt_condenseq_include'
require 'gt_consensus_sa_include'
require 'gt_csa_include'
require 'gt_csr_include.rb'
require 'gt_encseq_include'
require 'gt_eval_include'
require 'gt_extractfeat_include'
require 'gt_fastq_sample_include'
require 'gt_featureindex_include'
require 'gt_fingerprint_include'
require 'gt_genomediff_include'
require 'gt_gff3_include'
require 'gt_gff3validator_include'
require 'gt_gtf_to_gff3_include'
require 'gt_hop_include'
require 'gt_id_to_md5_include'
require 'gt_include'
require 'gt_inlineseq_include'
require 'gt_interfeat_include'
require 'gt_kmer_database_include'
require 'gt_linspace_align_include'
require 'gt_loccheck_include'
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
if python_tests_runnable? then
  require 'gt_python_include'
end
require 'gt_readjoiner_include'
require 'gt_readreads_include'
require 'gt_regioncov_include'
if ruby_tests_runnable? then
  require 'gt_ruby_include'
end
require 'gt_sambam_include'
require 'gt_script_filter_include'
require 'gt_scripts_include'
require 'gt_seed_extend_include'
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
require 'gt_tirvish_include'
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
require 'gt_repfind_include'
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

#output seed
puts "seed=#{$SEED}"

#start tests
$testsuite.run
