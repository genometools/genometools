#!/usr/bin/env ruby
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

#require 'escape'
require 'tmpdir'

class LTRharvestRun
  attr_accessor :config, :opts

  def initialize(job, config=nil, opts=nil)
    #default options, can be overridden using the 'config' parameter
    @config = {
      :seed =>"40",          #e.g. "40"
      :minlenltr => 100,     #e.g. "100"
      :maxlenltr => 1000,    #e.g. "1000"
      :mindistltr => 1500,   #e.g. "1500"
      :maxdistltr => 15000,  #e.g. "15000"
      :similar => 90,        #e.g. "90.0"
      :mintsd => 4,          #e.g. "4"
      :maxtsd => 20,         #e.g. "20"
      :motif => "",          #e.g. "tgca" ; not set by default
      :motifmis => "",       #e.g. "0" ; requires motif
      :overlaps => "best",   #e.g. "no"|"best"|"all"
      :xdrop => 5,           #e.g. "5"
      :mat =>  2,            #e.g. "2"
      :mis => -2,            #e.g. "-2"
      :ins => -3,            #e.g. "-3"
      :del => -3,            #e.g. "-3"
    }
    #default options, can be overridden using the 'opts' parameter
    #e.g. when files should be written to persistent directories instead of
    #tempdir
    #at least specify a path to the 'gt' executable!
    @opts ={
      :gtpath => nil,
      :outdir => Dir.tmpdir,
      :with_innerfile=>true,
      :with_outfile=>true,
      :with_gff3file=>true,
      :with_fastafile=>true
    }
    #job id, e.g. species
    @job = job
    #override defaults
    @config.merge!(config) unless config.nil?
    @opts.merge!(opts) unless opts.nil?
    @ids = []
    @seqs = {}
  end

  def add_seq(chr, file, config={})
    raise "sequence with id #{chr} duplicate" if @seqs.has_key?(chr)
    raise "file #{file}.prj not accessible" if !File.exist?("#{file}.prj")
    @ids.push(chr)
    @seqs[chr] = {
                   :chr => chr, :file => file, :config => config,
                   :innerfile => "#{@opts[:outdir]}/#{@job}_#{chr}"\
                                +"-innerregion.fsa",
                   :gff3file => "#{@opts[:outdir]}/#{@job}_#{chr}.gff3",
                   :fastafile => "#{@opts[:outdir]}/#{@job}_#{chr}.fsa",
                   :resultfile => "#{@opts[:outdir]}/#{@job}_#{chr}.result"
                 }
  end

  def each_seq
    @ids.each do |id|
      nr = id
      chr = @seqs[id]
      args = []
      if File.exist?("#{chr[:file]}.prj")
        args.push("ltrharvest")
        myconf = @config.merge(chr[:config]).reject{|k,v| (v.nil?||v=="")}
        myconf.each do |k,v|
          args.push("-#{k} #{v.to_s}")
        end
        args.push("-v")
        args.push("-out #{chr[:fastafile]}") unless !@opts[:with_fastafile]
        args.push("-gff3 #{chr[:gff3file]}") unless !@opts[:with_gff3file]
        args.push("-outinner #{chr[:innerfile]}") unless !@opts[:with_innerfile]
        args.push("-index #{chr[:file]}")
        args.push("> #{chr[:resultfile]}") unless !@opts[:with_outfile]
        yield nr, chr, args, myconf
      end
    end
  end

  # set 'test' to true to enable stest return code checking
  def run_seq(test=false)
    if @opts[:gtpath].nil?\
       || !File.exist?(@opts[:gtpath])\
       || !File.executable?(@opts[:gtpath]) then
      raise "gt binary not found or executable: #{@opts[:gtpath]}"
    end
    each_seq do |nr, chr, arglist, chr_cfg|
      if test then
        run "#{@opts[:gtpath]} #{arglist.join(" ")}", :maxtime => 500
      else
        STDERR.puts "Running #{@job}: seq '#{nr}'"
        success = Kernel.system("#{@opts[:gtpath]} #{arglist.join(" ")}")
      end
      if success or test then
        yield nr, chr[:resultfile], chr[:innerfile], \
              chr[:gff3file], chr[:fastafile]
      else
        raise "canceled command: #{@opts[:gtpath]} #{arglist.join(" ")}"
      end
    end
  end
end
