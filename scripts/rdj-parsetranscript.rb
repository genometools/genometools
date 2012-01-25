#! /usr/bin/env ruby
#
# Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
# Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg
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

## Extensions to base classes ##

class Float
  # round a float to a given n of decimals
  # or an integer if precision is 0
  def with_precision(precision)
    if precision > 0
      return ((self * (10 ** precision)).round.to_f) / (10 ** precision)
    else
      return self.round
    end
  end
end

class Fixnum
  # given a number of seconds, return a string with hours and minutes
  def s2hm
    h = self / 3600
    m = (self % 3600) / 60
    (h > 0) ? "#{h} h #{m} min" : "#{m} min"
  end
end

class String
  def escape_for_latex!
    self.gsub!("%", "\\%")
    self.gsub!("_", "\\_")
    return self
  end
end

# modified from: http://blog.sethladd.com/2007/03/creating-combinations-of.html
class Array
  def choose(k)
    return [[]] if nil? || empty? && k == 0
    return [] if nil? || empty? && k > 0
    return [[]] if size > 0 && k == 0
    c2 = clone
    c2.pop
    new_element = clone.pop
    return c2.choose(k) + c2.choose(k-1).map {|l| l << new_element}
  end
end

# for the -raw option:
# modified from: http://snippets.dzone.com/posts/show/5811
require 'yaml'
class Hash
  def to_yaml( opts = {} )
    YAML::quick_emit( object_id, opts ) do |out|
      out.map( taguri, to_yaml_style ) do |map|
        sort_by{|k, v| k.to_s}.each do |k, v|
          map.add( k, v )
        end
      end
    end
  end
end

class Numeric
  MemUnits = [:b, :kb, :Mb, :Gb, :Tb]
  MemUnitsDefaultPrecision = {:b => 0, :kb => 0, :Mb => 0, :Gb => 1, :Tb => 3}
  def mem_convert(from, to, precision = nil)
    if !MemUnits.include?(from) or !MemUnits.include?(to)
      raise ArgumentError, "allowed values: #{MemUnits.join(', ')}"
    end
    if precision.nil?
      precision = MemUnitsDefaultPrecision[to]
    end
    from_factor = MemUnits.index(from)
    to_factor = MemUnits.index(to)
    result = self.to_f * (1024**(from_factor)) / (1024**(to_factor))
    return result.with_precision(precision)
  end
end

## OutputParser: search information in different kinds of logs ##
module OutputParser

  module Common

    # output a warning if the given variable was
    # already assigned to a different value
    def warn_if_not_nil(var, newvalue)
      oldvalue = instance_variable_get(var)
      if !oldvalue.nil? && oldvalue != newvalue && !@warned
        @warned = true
        STDERR.puts "#{self.class} warning: multiple measurements found! "+
          "The value of #{var} was #{oldvalue.inspect}."
      end
    end
    private :warn_if_not_nil

    def disable_float_conversion
      @noautofloat = true
    end
    private :disable_float_conversion

    def enable_float_conversion
      @noautofloat = false
    end
    private :enable_float_conversion

    # assign one or more variables using a regex if a match is found
    def detect(line, regex, *values)
      detected = false
      @operation ||= []
      if line =~ regex
        detected = true
        values.each_with_index do |var, i|
          @operation[i+1] ||= ".to_f" unless @noautofloat
          value = eval("$#{i+1}#{@operation[i+1]}")
          warn_if_not_nil(var, value)
          instance_variable_set(var, value)
        end
      end
      @operation = []
      return detected
    end
    private :detect

    # autodiscover all public methods and return an hash with the values
    def results
      mets = self.class.instance_methods - ["results"] -
             self.class.superclass.instance_methods
      retval = {}
      mets.each do |met|
        retval[met.to_sym] = send(met)
      end
      return retval
    end

    def clean_readset_name
      if @readset
        @readset = File.basename(@readset)
        if @readset =~ /^([\w-]+).reads.fas$/
          @readset = $1
        end
      end
    end

    # --- class methods ---

    def self.included(klass)
      klass.extend ClassMethods
    end

    module ClassMethods

      # delegate a method to an instance variable
      def delegate(method, var, var_method = method)
        class_eval <<-enddef
          def #{method}
            @#{var}.#{var_method}
          end
        enddef
      end

    end

  end

  ### === Log parts === ###

  class Info
    include Common

    def initialize(log)
      @operation = []
      parse(log)
    end

    attr_reader :hostname, :compiler, :lastcommit

    private

    def parse(log)
      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /^# hostname: (.*)\s*$/, "@hostname")
      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /^# compiler: (\w+).*$/, "@compiler")
      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /^# last commit: (\w{7}).*$/, "@lastcommit")
    end
  end

  # the "parser" classes return time in s and space in Mb, as numerics;
  # lengths and counts are return as numeric values

  # reads most common time output formats (built-in time, GNU time, ...);
  # output a warning if more than one measurement is found in a log;
  # the *_time methods return nil or the time in seconds as float;
  # the cpu_time is defined as sum of user and sys;
  # cpu/sys/user is checked against real for plausibility (for n threads)
  class Time
    include Common

    def initialize(log, threads = 1)
      parse(log)
      @threads = threads
    end

    attr_reader :real_time

    def user_time
      if @user_time
        warn_if_too_large("user_time")
        return @user_time
      end
    end

    def sys_time
      if @sys_time
        warn_if_too_large("sys_time")
       return @sys_time
      end
    end

    def cpu_time
      if @user_time.nil? || @sys_time.nil?
        return nil
      else
        @cpu_time = @user_time.to_f + @sys_time.to_f
        warn_if_too_large("cpu_time")
        return @cpu_time
      end
    end

    def real_cpu_time_diff
      if @real_time.nil? || cpu_time.nil?
        return nil
      else
        return @real_time - cpu_time
      end
    end

    private

    Tolerance = 2 # seconds

    def warn_if_too_large(var)
      value = eval("@#{var}")
      if (!value.nil?) && (!@real_time.nil?)
        if value > (@real_time * @threads + Tolerance) && !@time_warning
          thrstr = ''
          thrstr = " * #{@threads}" if @threads > 1
          STDERR.puts "#{self.class} warning: #{var} > real#{thrstr}!"
          @time_warning = true
        end
      end
    end

    def detect_hms(line, regex, var)
      @t1, @t2, @t3 = nil, nil, nil
      detected = detect(line, regex, "@t1", "@t2", "@t3")
      if detected
        t = @t1 * 3600 + @t2 * 60 + @t3
        warn_if_not_nil(var, t)
        eval("#{var} = t")
      end
      return detected
    end

    def detect_ms(line, regex, var)
      @t1, @t2 = nil, nil
      detected = detect(line, regex, "@t1", "@t2")
      if detected
        t = @t1 * 60 + @t2
        warn_if_not_nil(var, t)
        eval("#{var} = t")
      end
      return detected
    end

    def detect_usr(line, us, r_with_h, r_without_h)
      detected = detect(line, us, "@user_time", "@sys_time")
      if detected
        detected2 = detect_hms(line, / (\d+):(\d+):(\d+.?\d*)elapsed/,
                              "@real_time")
        if !detected2
          detect_ms(line, / (\d+):(\d+.?\d*)elapsed/, "@real_time")
        end
      end
      return detected
    end

    def parse(log)
      log.each_line do |line|
        detect(line, /^user (\d+\.\d+)/, "@user_time")
        detect_ms(line, /^user\s+(\d+)m(\d+.\d+)s/, "@user_time")
        detect(line, /^real (\d+\.\d+)/, "@real_time" )
        detect_ms(line, /^real\s+(\d+)m(\d+.\d+)s/, "@real_time")
        detect(line, /^sys (\d+\.\d+)/, "@sys_time")
        detect(line, /^sys\s+(\d+)m(\d+.\d+)s/, "@sys_time")
        detect_usr(line, /^(\d+\.\d+)user (\d+\.\d+)system /,
                   / (\d+):(\d+):(\d+.?\d*)elapsed/,
                   / (\d+):(\d+.?\d*)elapsed/)
        detect_usr(line, /^(\d+\.\d+)u (\d+\.\d+)s/, / (\d+):(\d+):(\d+.\d+)/,
                   / (\d+):(\d+\.\d+)/)
      end
    end

  end

  # reads the output format of the measure_space_peak script
  # the values are returned as Floats and are expressed in Mb
  class MeasureSpacePeak
    include Common

    def initialize(log)
      parse(log)
      correct_units
    end

    attr_reader :vm_peak, :vm_size, :vm_lck, :vm_hwm, :vm_rss, :vm_data,
                :vm_stk, :vm_exe, :vm_lib, :vm_pte

    private

    def parse(log)
      log.each_line do |line|
        detect(line, /VmPeak:\s+(\d+) kB/, "@vm_peak")
        detect(line, /VmSize:\s+(\d+) kB/, "@vm_size")
        detect(line, /VmLck:\s+(\d+) kB/, "@vm_lck")
        detect(line, /VmHWM:\s+(\d+) kB/, "@vm_hwm")
        detect(line, /VmRSS:\s+(\d+) kB/, "@vm_rss")
        detect(line, /VmData:\s+(\d+) kB/, "@vm_data")
        detect(line, /VmStk:\s+(\d+) kB/, "@vm_stk")
        detect(line, /VmExe:\s+(\d+) kB/, "@vm_exe")
        detect(line, /VmLib:\s+(\d+) kB/, "@vm_lib")
        detect(line, /VmPTE:\s+(\d+) kB/, "@vm_pte")
      end
    end

    def correct_units
      @vm_peak = @vm_peak.mem_convert(:kb, :Mb) if @vm_peak
      @vm_size = @vm_size.mem_convert(:kb, :Mb) if @vm_size
      @vm_lck  = @vm_lck.mem_convert(:kb, :Mb)  if @vm_lck
      @vm_hwm  = @vm_hwm.mem_convert(:kb, :Mb)  if @vm_hwm
      @vm_rss  = @vm_rss.mem_convert(:kb, :Mb)  if @vm_rss
      @vm_data = @vm_data.mem_convert(:kb, :Mb) if @vm_data
      @vm_stk  = @vm_stk.mem_convert(:kb, :Mb)  if @vm_stk
      @vm_exe  = @vm_exe.mem_convert(:kb, :Mb)  if @vm_exe
      @vm_lib  = @vm_lib.mem_convert(:kb, :Mb)  if @vm_lib
      @vm_pte  = @vm_pte.mem_convert(:kb, :Mb)  if @vm_pte
    end

  end

  # parser for the genometools -showtime and -spacepeak output
  class Genometools
    include Common

    def initialize(log)
      @time = {:real => {}, :user => {}, :sys => {}}
      @space_peak = {}
      parse(log)
    end

    attr_reader :mmap_peak, :alloc_peak, :time, :space_peak

    def real_time_overall
      @time[:real]["overall"]
    end

    def user_time_overall
      @time[:user]["overall"]
    end

    def sys_time_overall
      @time[:sys]["overall"]
    end

    def cpu_time_overall
      if sys_time_overall && user_time_overall
        sys_time_overall + user_time_overall
      end
    end

    def real_cpu_time_diff_overall
      if real_time_overall.nil? || cpu_time_overall.nil?
        return nil
      else
        return real_time_overall - cpu_time_overall
      end
    end

    def space_peak_overall
      if @overall_peak
        return @overall_peak
      elsif @space_peak["overall"]
        return @space_peak["overall"]
      # otherwise use sum of mmap and alloc peaks
      elsif mmap_peak.nil? or alloc_peak.nil?
        return nil
      else
        return mmap_peak + alloc_peak
      end
    end

    private

    def parse(log)
      log.each_line do |line|
        detect(line, /# mmap space peak in megabytes: (\d+.\d+)/,
               "@mmap_peak")
        detect(line, /# space peak in megabytes: (\d+.\d+)/,
               "@alloc_peak")
        detect(line, /# combined space peak in megabytes: (\d+.\d+)/,
               "@overall_peak")
        @operation[1] = "" #next detect: no auto-conversion into Float for $1
        d = detect(line, /# TIME (.*?) (\d+\.\d+)/, "@tmpdesc", "@tmptime")
        if d
          @time[:real][@tmpdesc] = @tmptime
          d = detect(line, /\(user: (\d+\.\d+); sys: (\d+\.\d+)\)/, "@tmputime",
                     "@tmpstime")
          if d
            @time[:user][@tmpdesc] = @tmputime
            @time[:sys][@tmpdesc] = @tmpstime
          end
          @tmpdesc, @tmptime, @tmputime, @tmpstime = nil, nil, nil, nil
        end
        @operation[1] = "" #next detect: no auto-conversion into Float for $1
        d = detect(line, /# SPACE PEAK (.*?) (\d+\.\d+) Mb/, "@tmpdesc",
                   "@tmpspace")
        if d
          @space_peak[@tmpdesc] = @tmpspace
          @tmpdesc, @tmpspace = nil, nil
        end
      end
    end

  end

  class SeqstatContigs
    include Common

    def initialize(log)
      parse(log)
    end

    attr_reader :n50, :n80, :lngst, :nofcontigs, :tcontigslen

    private

    def parse(log)
      detect(log, /# N50:\s+(\d+)/, "@n50")
      detect(log, /# N80:\s+(\d+)/, "@n80")
      detect(log, /# longest contig:\s+(\d+)/, "@lngst")
      detect(log,
           /# number of contigs:\s+(\d+)\n# total length:\s+(\d+)/m,
           "@nofcontigs", "@tcontigslen")
      if !@tcontigslen
        detect(log,
             /# number of contigs:\s+(\d+)\n# total contigs length:\s+(\d+)/m,
             "@nofcontigs", "@tcontigslen")
      end
      if log =~ /# no contigs/
        @nofcontigs = 0
        @n50 = 0
        @n80 = 0
        @lngst = 0
        @tcontigslen = 0
      end
    end
  end

  ### ===== programs ===== ###

  class SgaIndex
    include Common

    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
      parse(log)
    end

    attr_reader :sga_index_d, :nofreads, :readset

    delegate :index_hostname,      :i,   :hostname
    delegate :index_real_time,     :t,   :real_time
    delegate :index_user_time,     :t,   :user_time
    delegate :index_sys_time,      :t,   :sys_time
    delegate :index_cpu_time,      :t,   :cpu_time
    delegate :index_real_cpu_time_diff, :t,   :real_cpu_time_diff
    delegate :index_hwm_space,     :msp, :vm_hwm
    delegate :index_vm_peak,       :msp, :vm_peak

    def index_pipeline
      "Sga"
    end

    private

    def parse(log)
      detect(log, /sga index.*-d (\d+).*/, "@sga_index_d")
      arr = log.scan(/Merge\d+: \[\d+,(\d+)\].*/)
      if !arr.empty?
        @nofreads = arr.flatten.map{|x|x.to_i}.max
      end
      disable_float_conversion
      detect(log, /Building index for (\S+)/, "@readset")
      clean_readset_name
    end

  end

  class SgaOverlap
    include Common

    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
      parse(log)
    end

    attr_reader :minlen

    delegate :overlap_hostname,      :i,   :hostname
    delegate :overlap_real_time,     :t,   :real_time
    delegate :overlap_user_time,     :t,   :user_time
    delegate :overlap_sys_time,      :t,   :sys_time
    delegate :overlap_cpu_time,      :t,   :cpu_time
    delegate :overlap_real_cpu_time_diff, :t,   :real_cpu_time_diff
    delegate :overlap_hwm_space,     :msp, :vm_hwm
    delegate :overlap_vm_peak,       :msp, :vm_peak

    def overlap_pipeline
      "Sga"
    end

    private

    def parse(log)
      detect(log, /.*sga overlap.*?-m (\d+).*/, "@minlen")
    end

  end

  class SgaAssembly
    include Common

    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
      @ssc = ::OutputParser::SeqstatContigs.new(log)
    end

    delegate :assembly_hostname,      :i,   :hostname
    delegate :assembly_real_time,     :t,   :real_time
    delegate :assembly_user_time,     :t,   :user_time
    delegate :assembly_sys_time,      :t,   :sys_time
    delegate :assembly_cpu_time,      :t,   :cpu_time
    delegate :assembly_real_cpu_time_diff, :t,   :real_cpu_time_diff
    delegate :assembly_hwm_space,     :msp, :vm_hwm
    delegate :assembly_vm_peak,       :msp, :vm_peak
    delegate :nofcontigs,             :ssc, :nofcontigs
    delegate :lngst,                  :ssc, :lngst
    delegate :n50,                    :ssc, :n50
    delegate :n80,                    :ssc, :n80
    delegate :tcontigslen,            :ssc, :tcontigslen

    def assembly_pipeline
      "Sga"
    end

  end

  class Leap
    include Common

    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
      @ssc = ::OutputParser::SeqstatContigs.new(log)
      parse(log)
    end

    attr_reader :minlen, :tcontigslen, :lngst, :nofcontigs, :nofreads,
      :contreads, :contreads_percent, :tlen, :readset

    delegate :overall_hostname,      :i,   :hostname
    delegate :overall_real_time,     :t,   :real_time
    delegate :overall_user_time,     :t,   :user_time
    delegate :overall_sys_time,      :t,   :sys_time
    delegate :overall_cpu_time,      :t,   :cpu_time
    delegate :overall_real_cpu_time_diff, :t,   :real_cpu_time_diff
    delegate :overall_hwm_space,     :msp, :vm_hwm
    delegate :overall_vm_peak,       :msp, :vm_peak
    delegate :n50,                   :ssc, :n50
    delegate :n80,                   :ssc, :n80

    def overlap_pipeline
      "Leap"
    end

    def readlen
      if @minreadlen == @maxreadlen
        @minreadlen
      else
        "variable"
      end
    end

    private

    def parse(log)
      detect(log, /^overlap threshold: (\d+)/, "@minlen")
      detect(log, /^total contig length: (\d+)/, "@tcontigslen")
      detect(log, /^max contig length: (\d+)/, "@lngst")
      detect(log, /^#contigs: (\d+)/, "@nofcontigs")
      detect(log, /^#reads = (\d+)/, "@nofreads")
      detect(log, /remove contained reads:.*\n#reads = (\d+)/, "@nofcfreads")
      if !@nofcfreads.nil?
        @contreads = @nofreads - @nofcfreads
        @contreads_percent = @contreads * 100 / @nofreads
      end
      detect(log, /^#bases = (\d+)/, "@tlen")
      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /^read input file (\S+)/, "@readset")
      clean_readset_name
      detect(log, /^minReadLen = (\d+)/, "@minreadlen")
      detect(log, /^maxReadLen = (\d+)/, "@maxreadlen")
    end
  end

  module EdenaCommon

    def contreads
      return nil if @notamb.nil? or @nofnodes.nil?
      return @notamb.to_i - @nofnodes.to_i
    end

    def contreads_percent
      return nil if contreads.nil?
      return contreads.to_f / @nofreads.to_f * 100
    end

  end

  class EdenaOverlap
    include Common
    include EdenaCommon

    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
      parse(log)
    end

    attr_reader :minlen, :nofreads, :introvl, :readlen, :readset

    def tlen
      if (readlen && nofreads)
        return readlen * nofreads
      end
    end

    delegate :overlap_hostname,      :i,   :hostname
    delegate :overlap_real_time,     :t,   :real_time
    delegate :overlap_user_time,     :t,   :user_time
    delegate :overlap_sys_time,      :t,   :sys_time
    delegate :overlap_cpu_time,      :t,   :cpu_time
    delegate :overlap_real_cpu_time_diff, :t,   :real_cpu_time_diff
    delegate :overlap_hwm_space,     :msp, :vm_hwm
    delegate :overlap_vm_peak,       :msp, :vm_peak

    def overlap_pipeline
      "Edena"
    end

    private

    def parse(log)
      log.each_line do |line|
        if @nextisreadset
          detect(line, /(\S+)/, "@readset")
          @nextisreadset = false
        end
        detect(line, /Opening file (\S+)/, "@readset")
        detect(line, /Opening file/, "@nextisreadset") if @readset.nil?
        detect(line, /reads length:\s+(\d+)/i, "@readlen")
        detect(line, /number of nodes:\s+(\d+)/, "@nofnodes")
        detect(line, /number of edges:\s+(\d+)/, "@introvl")
        detect(line, /minimum overlap size:\s+(\d+)/, "@minlen")
        detect(line, /fasta entries:\s+(\d+)/, "@nofreads")
        detect(line, /sequences OK:\s+(\d+)/, "@notamb")
        detect(line, /Number of reads:\s+(\d+)/, "@nofreads")
        @operation[1] = "" #next detect: no auto-conversion into Float for $1
        clean_readset_name
      end
     end
  end

  class EdenaAssembly
    include Common
    include EdenaCommon

    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
      @ssc = ::OutputParser::SeqstatContigs.new(log)
      parse(log)
    end

    attr_reader :minlen, :lngst, :n50, :nofreads, :introvl, :nofcontigs,
      :tcontigslen

    delegate :assembly_hostname,      :i,   :hostname
    delegate :assembly_real_time,     :t,   :real_time
    delegate :assembly_user_time,     :t,   :user_time
    delegate :assembly_sys_time,      :t,   :sys_time
    delegate :assembly_cpu_time,      :t,   :cpu_time
    delegate :assembly_real_cpu_time_diff, :t,   :real_cpu_time_diff
    delegate :assembly_hwm_space,     :msp, :vm_hwm
    delegate :assembly_vm_peak,       :msp, :vm_peak
    delegate :ssc_nofcontigs,         :ssc, :nofcontigs
    delegate :n80,                    :ssc, :n80
    delegate :ssc_tcontigslen,        :ssc, :tcontigslen

    def assembly_pipeline
      "Edena"
    end

    private

    def parse(log)
      log.each_line do |line|
        detect(line, /^\s*N50:\s+(\d+)/, "@n50")
        detect(line, /max:\s+(\d+)/, "@lngst")
        detect(line, /minimum overlap size:\s+(\d+)/, "@minlen")
        detect(line, /number of nodes:\s+(\d+)/, "@nofnodes") if @nofnodes.nil?
        detect(line, /number of edges:\s+(\d+)/, "@introvl")
        detect(line, /number of reads:\s+(\d+)/, "@nofreads")
        detect(line, /Number of contigs:\s+(\d+)/, "@nofcontigs")
        @nofcontigs ||= ssc_nofcontigs
        detect(line, /Total number of bases:\s+(\d+)/, "@tcontigslen")
        @tcontigslen ||= ssc_tcontigslen
      end
    end
  end

  class Prefilter
    include Common
    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @gt  = ::OutputParser::Genometools.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
      parse(log)
    end

    attr_reader :prefilter_indexname, :nofreads, :tlen, :contreads,
      :contreads_percent, :readset, :readlen

    delegate :prefilter_hostname,      :i,   :hostname
    delegate :prefilter_compiler,      :i,   :compiler
    delegate :prefilter_lastcommit,    :i,   :lastcommit
    delegate :prefilter_gt_real_time,  :gt,  :real_time_overall
    delegate :prefilter_gt_user_time,  :gt,  :user_time_overall
    delegate :prefilter_gt_sys_time,   :gt,  :sys_time_overall
    delegate :prefilter_gt_cpu_time,   :gt,  :cpu_time_overall
    delegate :prefilter_real_cpu_time_diff, :gt,  :real_cpu_time_diff_overall
    delegate :prefilter_gt_space,      :gt,  :space_peak_overall
    delegate :prefilter_alloc,         :gt,  :alloc_peak
    delegate :prefilter_mmap,          :gt,  :mmap_peak
    delegate :prefilter_real_time,     :t,   :real_time
    delegate :prefilter_user_time,     :t,   :user_time
    delegate :prefilter_sys_time,      :t,   :sys_time
    delegate :prefilter_cpu_time,      :t,   :cpu_time
    delegate :prefilter_real_cpu_time_diff, :t,   :real_cpu_time_diff
    delegate :prefilter_hwm_space,     :msp, :vm_hwm
    delegate :prefilter_vm_peak,       :msp, :vm_peak

    def prefilter_pipeline
      "Rdj"
    end

    private

    def parse(log)
      @operation = []
      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /readset name = (\S+)/, "@readset")
      clean_readset_name
      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /indexname = (\S+)/, "@prefilter_indexname")
      detect(log, /total input length = (\d+)/, "@tlen")
      detect(log, /total length of complete readset = (\d+)/, "@tlen")
      detect(log, /read length = (\d+)/, "@readlen")
      detect(log, /input sequences = (\d+)/, "@nofreads")
      detect(log, /number of reads in complete readset = (\d+)/, "@nofreads")
      detect(log, /contained sequences = (\d+)/, "@contreads")
      detect(log, /contained reads = (\d+)/, "@contreads")
      if @contreads && @nofreads
        @contreads_percent = @contreads * 100 / @nofreads
      end
    end
  end

  class Sequniq
    include Common
    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @gt  = ::OutputParser::Genometools.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
    end
    delegate :prefilter_hostname,      :i,   :hostname
    delegate :prefilter_compiler,      :i,   :compiler
    delegate :prefilter_lastcommit,    :i,   :lastcommit
    delegate :prefilter_gt_real_time,  :gt,  :real_time_overall
    delegate :prefilter_gt_user_time,  :gt,  :user_time_overall
    delegate :prefilter_gt_sys_time,   :gt,  :sys_time_overall
    delegate :prefilter_gt_cpu_time,   :gt,  :cpu_time_overall
    delegate :prefilter_real_cpu_time_diff, :gt,  :real_cpu_time_diff_overall
    delegate :prefilter_gt_space,      :gt,  :space_peak_overall
    delegate :prefilter_alloc,         :gt,  :alloc_peak
    delegate :prefilter_mmap,          :gt,  :mmap_peak
    delegate :prefilter_real_time,     :t,   :real_time
    delegate :prefilter_user_time,     :t,   :user_time
    delegate :prefilter_sys_time,      :t,   :sys_time
    delegate :prefilter_cpu_time,      :t,   :cpu_time
    delegate :prefilter_real_cpu_time_diff, :t,   :real_cpu_time_diff
    delegate :prefilter_hwm_space,     :msp, :vm_hwm
    delegate :prefilter_vm_peak,       :msp, :vm_peak

    def prefilter_pipeline
      "Rdj"
    end

  end

  class Suffixerator
    include Common

    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @gt  = ::OutputParser::Genometools.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
      parse(log)
    end

    attr_reader :parts, :prefixlength, :maxinsertionsort, :maxbltriesort,
                :maxcountingsort, :relevant, :suftabuint, :memlimit,
                :sat, :sfxtotallength, :index_indexname, :index_bits

    delegate :index_hostname,      :i,   :hostname
    delegate :index_compiler,      :i,   :compiler
    delegate :index_lastcommit,    :i,   :lastcommit
    delegate :index_gt_real_time,  :gt,  :real_time_overall
    delegate :index_gt_user_time,  :gt,  :user_time_overall
    delegate :index_gt_sys_time,   :gt,  :sys_time_overall
    delegate :index_gt_cpu_time,   :gt,  :cpu_time_overall
    delegate :index_real_cpu_time_diff, :gt,  :real_cpu_time_diff_overall
    delegate :index_gt_space,      :gt,  :space_peak_overall
    delegate :index_alloc,         :gt,  :alloc_peak
    delegate :index_mmap,          :gt,  :mmap_peak
    delegate :index_real_time,     :t,   :real_time
    delegate :index_user_time,     :t,   :user_time
    delegate :index_sys_time,      :t,   :sys_time
    delegate :index_cpu_time,      :t,   :cpu_time
    delegate :index_real_cpu_time_diff, :t,   :real_cpu_time_diff
    delegate :index_hwm_space,     :msp, :vm_hwm
    delegate :index_vm_peak,       :msp, :vm_peak

    def index_hwm_gtsp_delta
      return nil if !index_hwm_space or !index_gt_space
      return (index_hwm_space - index_gt_space) * 100 / index_gt_space
    end

    def index_pipeline
      "Rdj"
    end

    private

    def parse(log)
      detect(log, /-parts (\d+)/, "@parts")
      detect(log, /# parts=(\d+)/, "@parts")
      detect(log, /# derived parts=(\d+)/, "@parts")
      detect(log, /# prefixlength=(\d+)/, "@prefixlength")
      if @prefixlength.nil?
        detect(log, /# automatically determined prefixlength=(\d+)/,
               "@prefixlength")
      end
      detect(log, /# maxinsertionsort=(\d+)/, "@maxinsertionsort")
      detect(log, /# maxbltriesort=(\d+)/, "@maxbltriesort")
      detect(log, /# maxcountingsort=(\d+)/, "@maxcountingsort")
      detect(log, /# relevant suffixes=(\d+.\d+)\%/, "@relevant")
      if @relevant.nil?
        @relevant = "spmopt off"
      end
      @suftabuint = (log =~ /suftab uses 32bit values/) ? "on" : "off"
      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /# maximumspace=(.*)\s*\n/, "@memlimit")
      detect(log, /# totallength=(\d+)/, "@sfxtotallength")
      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /# indexname=\"(.*)\"\s*\n/, "@index_indexname")
      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /# sat=(.*)\s*\n/, "@sat")
      if @sat.nil?
        @operation[1] = "" #next detect: no auto-conversion into Float for $1
        detect(log, /# init character encoding \((\S+)/, "@sat")
      end
      detect(log, /# sizeof \(unsigned long\)=(\d+)/, "@index_bits")
    end

  end

  class RdjOverlap
    include Common

    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @gt  = ::OutputParser::Genometools.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
      @index = ::OutputParser::Suffixerator.new(log)
      parse(log)
    end

    attr_reader :minlen, :trovl, :introvl, :contreads, :contreads_percent,
                :parts, :nofreads, :overlap_indexname

    delegate :overlap_hostname,         :i,   :hostname
    delegate :overlap_compiler,         :i,   :compiler
    delegate :overlap_lastcommit,       :i,   :lastcommit
    delegate :overlap_gt_real_time,     :gt,  :real_time_overall
    delegate :overlap_gt_user_time,     :gt,  :user_time_overall
    delegate :overlap_gt_sys_time,      :gt,  :sys_time_overall
    delegate :overlap_gt_cpu_time,      :gt,  :cpu_time_overall
    delegate :overlap_gt_real_cpu_time_diff, :gt,  :real_cpu_time_diff_overall
    delegate :overlap_gt_space,         :gt,  :space_peak_overall
    delegate :overlap_alloc,            :gt,  :alloc_peak
    delegate :overlap_mmap,             :gt,  :mmap_peak
    delegate :overlap_real_time,        :t,   :real_time
    delegate :overlap_user_time,        :t,   :user_time
    delegate :overlap_sys_time,         :t,   :sys_time
    delegate :overlap_cpu_time,         :t,   :cpu_time
    delegate :overlap_real_cpu_time_diff,    :t,   :real_cpu_time_diff
    delegate :overlap_hwm_space,        :msp, :vm_hwm
    delegate :overlap_vm_peak,          :msp, :vm_peak

    def overlap_hwm_gtsp_delta
      return nil if !overlap_hwm_space or !overlap_gt_space
      return (overlap_hwm_space - overlap_gt_space) * 100 / overlap_gt_space
    end

    def overlap_pipeline
      "Rdj"
    end

    private

    def parse(log)
      detect(log, /-parts (\d+)/, "@parts")
      detect(log, /derived parts=(\d+)/, "@parts")
      detect(log, /minimal spm length.*?(\d+)/, "@minlen")
      detect(log, /minimal match length.*?(\d+)/, "@minlen")
      if @minlen.nil?
        detect(log, /minimal overlap.*?(\d+)/, "@minlen")
      end
      detect(log, /omitted transitive suffix-prefix matches: (\d+)/, "@trovl")
      detect(log, /number of transitive suffix-prefix matches = (\d+)/,
             "@trovl")
      detect(log, /-- suffix-prefix matches: (\d+)/, "@introvl")
      detect(log, /number of irreducible suffix-prefix matches = (\d+)/,
             "@introvl")
      detect(log, /contained: (\d+) \((\d+.\d+) %\)/, "@contreads",
             "@contreads_percent")
      detect(log, /# reads: (\d+)/, "@nofreads")
      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /indexname = \"(\S+)\"/, "@overlap_indexname")
    end

  end

  class RdjAssembly
    include Common

    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @ssc = ::OutputParser::SeqstatContigs.new(log)
      @gt  = ::OutputParser::Genometools.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
      parse(log)
    end

    attr_reader :minlen, :nofreads, :assembly_indexname

    delegate :nofcontigs,                :ssc, :nofcontigs
    delegate :lngst,                     :ssc, :lngst
    delegate :n50,                       :ssc, :n50
    delegate :n80,                       :ssc, :n80
    delegate :tcontigslen,               :ssc, :tcontigslen
    delegate :assembly_gt_real_time,     :gt,  :real_time_overall
    delegate :assembly_gt_user_time,     :gt,  :user_time_overall
    delegate :assembly_gt_sys_time,      :gt,  :sys_time_overall
    delegate :assembly_gt_cpu_time,      :gt,  :cpu_time_overall
    delegate :assembly_gt_real_cpu_time_diff, :gt,  :real_cpu_time_diff_overall
    delegate :assembly_real_time,        :t,   :real_time
    delegate :assembly_user_time,        :t,   :user_time
    delegate :assembly_sys_time,         :t,   :sys_time
    delegate :assembly_cpu_time,         :t,   :cpu_time
    delegate :assembly_real_cpu_time_diff,    :t,   :real_cpu_time_diff
    delegate :assembly_gt_space,         :gt,  :space_peak_overall
    delegate :assembly_alloc,            :gt,  :alloc_peak
    delegate :assembly_mmap,             :gt,  :mmap_peak
    delegate :assembly_hwm_space,        :msp, :vm_hwm
    delegate :assembly_vm_peak,          :msp, :vm_peak
    delegate :assembly_hostname,         :i,   :hostname
    delegate :assembly_compiler,         :i,   :compiler
    delegate :assembly_lastcommit,       :i,   :lastcommit

    def assembly_hwm_gtsp_delta
      return nil if !assembly_hwm_space or !assembly_gt_space
      return (assembly_hwm_space - assembly_gt_space) * 100 / assembly_gt_space
    end

    def assembly_pipeline
      "Rdj"
    end

    private

    def parse(log)
      detect(log, /minimal overlap.*?(\d+)/, "@minlen")
      detect(log, /# reads: (\d+)/, "@nofreads")
      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /indexname = \"(\S+)\"/, "@assembly_indexname")
    end

  end

  class RdjBoth
    include Common

    def initialize(log)
      @i   = ::OutputParser::Info.new(log)
      @gt  = ::OutputParser::Genometools.new(log)
      @t   = ::OutputParser::Time.new(log)
      @msp = ::OutputParser::MeasureSpacePeak.new(log)
      parse(log)
    end

    attr_reader :minlen, :trovl, :introvl, :contreads,
                :contreads_percent, :nofreads,
                :overlap_gt_space,
                :overlap_gt_real_time,
                :overlap_gt_user_time,
                :overlap_gt_sys_time,
                :overlap_gt_cpu_time,
                :overlap_gt_real_cpu_time_diff,
                :overlap_and_assembly_indexname

    attr_reader :n50, :lngst, :nofcontigs, :tcontigslen,
                :assembly_gt_space,
                :assembly_gt_real_time,
                :assembly_gt_user_time,
                :assembly_gt_sys_time,
                :assembly_gt_cpu_time,
                :assembly_gt_real_cpu_time_diff

    delegate :overlap_and_assembly_hostname,         :i,   :hostname
    delegate :overlap_and_assembly_compiler,         :i,   :compiler
    delegate :overlap_and_assembly_lastcommit,       :i,   :lastcommit
    delegate :overlap_and_assembly_gt_real_time,     :gt,  :real_time_overall
    delegate :overlap_and_assembly_gt_user_time,     :gt,  :user_time_overall
    delegate :overlap_and_assembly_gt_sys_time,      :gt,  :sys_time_overall
    delegate :overlap_and_assembly_gt_cpu_time,      :gt,  :cpu_time_overall
    delegate :overlap_and_assembly_gt_real_cpu_time_diff, :gt,
      :real_cpu_time_diff_overall
    delegate :overlap_and_assembly_gt_space,         :gt,  :space_peak_overall
    delegate :overlap_and_assembly_alloc,            :gt,  :alloc_peak
    delegate :overlap_and_assembly_mmap,             :gt,  :mmap_peak
    delegate :overlap_and_assembly_real_time,        :t,   :real_time
    delegate :overlap_and_assembly_user_time,        :t,   :user_time
    delegate :overlap_and_assembly_sys_time,         :t,   :sys_time
    delegate :overlap_and_assembly_cpu_time,         :t,   :cpu_time
    delegate :overlap_and_assembly_real_cpu_time_diff,    :t,
      :real_cpu_time_diff
    delegate :overlap_and_assembly_hwm_space,        :msp, :vm_hwm
    delegate :overlap_and_assembly_vm_peak,          :msp, :vm_peak

    def overlap_and_assembly_hwm_gtsp_delta
      if !overlap_and_assembly_hwm_space or !overlap_and_assembly_gt_space
        return nil
      end
      return (overlap_and_assembly_hwm_space - overlap_and_assembly_gt_space) *
        100 / overlap_and_assembly_gt_space
    end

    def overlap_and_assembly_pipeline
      "Rdj"
    end

    private

    StageDesc =
      {
       # Overlap phase
         :overlap_phase => "Overlap Phase",
         # containments
           :CNT1 => "remove contained reads", # old
           :CNT2 => "identify contained reads",
           :CNT3 => "save contained reads list",
           :CNT4 => "load contained reads list",
         # spm
           :SPM1 => "count overlaps for each reads", # old
           :SPM2 => "determine suffix-prefix matches",
       # Assembly Phase
         :assembly_phase => "Assembly Phase",
         # create graph
           :SG1 => "load overlap counts",
           :SG2 => "create overlap graph", # old
           :SG3 => "create string graph",
           :SG4 => "save string graph",
           :SG5 => "export string graph in GraphViz dot format",
           :SG6 => "output graph information as spm list",
           :SG7 => "output string graph adjacence list",
           :SG8 => "load string graph",
         # graph manipulation
           :GM1 => "sort edges by length",
           :GM2 => "reduce self-match edges",
           :GM3 => "reduce withrc-match edges",
           :GM4 => "reduce transitive edges",
           :GM5 => "compact edges list",
           :GM6 => "reduce non-maximal edges",
           :GM7 => "remove dead-end paths",
           :GM8 => "remove p-bubbles",
         # output contigs
           :OC1 => "output contigs", # old
           :OC2 => "traverse graph and output contigs",
           :OC3 => "convert symbolic contigs into fasta",
           :OC4 => "output contig sequences from paths",
      }

    PhaseStages =
      {:overlap  => [:overlap_phase, :CNT1, :CNT2, :CNT3, :CNT4, :SPM1, :SPM2],
       :assembly => [:assembly_phase, :SG1, :SG2, :SG3, :SG4, :SG5, :SG6, :SG7,
                     :SG8, :GM1, :GM2, :GM3, :GM4, :GM5, :GM6, :GM7, :GM8, :OC1,
                     :OC2, :OC3]}

    def phase_time(variety, phase)
      time = 0.0
      PhaseStages[phase].each do |stage|
        t = @gt.time[variety][StageDesc[stage]]
        (time += t) if t
      end
      return time if time > 0.0 || variety == :sys
    end

    def phase_space(phase)
      space = []
      PhaseStages[phase].each do |stage|
        s = @gt.space_peak[StageDesc[stage]]
        (space.push(s)) if s
      end
      space.compact!
      return space.max if space.size > 0
    end

    def parse(log)
      detect(log, /# N50:\s+(\d+)/, "@n50")
      detect(log, /# longest contig:\s+(\d+)/, "@lngst")
      detect(log, /minimal spm length.*?(\d+)/, "@minlen")
      if @minlen.nil?
        detect(log, /minimal overlap.*?(\d+)/, "@minlen")
      end
      detect(log, /# number of contigs:\s(\d+)\n# total length:\s+(\d+)/m,
             "@nofcontigs", "@tcontigslen")
      d = detect(log, /nof tr edges: (\d+)/, "@trovl")
      @trovl = @trovl/2 if d
      detect(log, /omitted transitive suffix-prefix matches: (\d+)/,
             "@trovl")
      detect(log, /valid overlaps: (\d+)/, "@introvl")
      detect(log, /-- suffix-prefix matches: (\d+)/, "@introvl")
      detect(log, /contained: (\d+) \((\d+.\d+) %\)/, "@contreads",
             "@contreads_percent")
      detect(log, /# reads: (\d+)/, "@nofreads")

      @overlap_gt_real_time = phase_time(:real, :overlap)
      @overlap_gt_user_time = phase_time(:user, :overlap)
      @overlap_gt_sys_time = phase_time(:sys, :overlap)
      if @overlap_gt_user_time and @overlap_gt_sys_time
        @overlap_gt_cpu_time = @overlap_gt_user_time + @overlap_gt_sys_time
      end
      if @overlap_gt_cpu_time and @overlap_gt_real_time
        @overlap_gt_real_cpu_time_diff = @overlap_gt_real_time -
          @overlap_gt_cpu_time
      end
      @overlap_gt_space = phase_space(:overlap)
      @assembly_gt_real_time = phase_time(:real, :assembly)
      @assembly_gt_user_time = phase_time(:user, :assembly)
      @assembly_gt_sys_time = phase_time(:sys, :assembly)
      if @assembly_gt_user_time and @assembly_gt_sys_time
        @assembly_gt_cpu_time = @assembly_gt_user_time + @assembly_gt_sys_time
      end
      if @assembly_gt_cpu_time and @assembly_gt_real_time
        @assembly_gt_real_cpu_time_diff =
           @assembly_gt_real_time - @assembly_gt_cpu_time
      end
      @assembly_gt_space = phase_space(:assembly)

      @operation[1] = "" #next detect: no auto-conversion into Float for $1
      detect(log, /indexname = \"(\S+)\"/, "@overlap_and_assembly_indexname")
    end

  end

  class ResultsFile
    include Common

    def initialize(log)
      parse(log)
      correct_units
    end

    attr_reader :minlen, :nofreads, :trovl, :introvl, :n50, :lngst, :parts,
                :contreads, :contreads_percent, :index_time, :index_space,
                :overlap_time, :overlap_space, :assembly_time

    private

    def parse(log)
      log.each_line do |line|
        detect(line, /Min SPM length\s+(\d+)/, "@minlen")
        detect(line, /Number of reads \(M\)\s+(\d+\.\d+)/, "@nofreads")
        detect(line, /Contained reads \(k\)\s+(\d+)/     , "@contreads")
        detect(line, /Contained %\s+(\d+\.\d+)/          , "@contreads_percent")
        detect(line, /Transitive SPM \(M\)\s+(\d+\.\d+)/ , "@trovl")
        detect(line, /Irreducible SPM \(M\)\s+(\d+\.\d+)/, "@introvl")
        detect(line, /N50 contig \(kbp\)\s+(\d+\.\d+)/   , "@n50")
        detect(line, /Lngst contig \(kbp\)\s+(\d+\.\d+)/ , "@lngst")
        detect(line, /Index time \(s\)\s+(\d+)/          , "@index_time")
        detect(line, /Index space \(Gb\)\s+(\d+\.\d+)/   , "@index_space")
        detect(line, /Overlap time \(s\)\s+(\d+)/        , "@overlap_time")
        detect(line, /Assembly time \(s\)\s+(\d+)/       , "@assembly_time")
        detect(line, /Overlap space \(Gb\)\s+(\d+\.\d+)/ , "@overlap_space")
        detect(line, /Assembly space \(Gb\)\s+(\d+\.\d+)/, "@assembly_space")
      end
    end

    def correct_units
      @nofreads       = @nofreads  * 1_000_000                if @nofreads
      @trovl          = @trovl     * 1_000_000                if @trovl
      @introvl        = @introvl   * 1_000_000                if @introvl
      @contreads      = @contreads * 1_000                    if @contreads
      @n50            = @n50       * 1_000                    if @n50
      @lngst          = @lngst     * 1_000                    if @lngst
      @index_space    = @index_space.mem_convert(:Gb, :Mb)    if @index_space
      @assembly_space = @assembly_space.mem_convert(:Gb, :Mb) if @assembly_space
    end

  end

end

## Table: data formatter for output ##

class Table

  LabelsColSpacing = 3
  DataColSpacing = 5

  def initialize
    @format = :txt
    @tex_page = :completedocument
  end

  def data=(bidimensional_array)
    @table = bidimensional_array
    @nofcolumns = @table.size > 0 ? @table.first.size : 0
    @labels_col_width = determine_labels_col_width
    @data_col_witdh = determine_data_col_width
  end

  attr_writer :format, :tex_page, :nosideways

  def show
    retval = ""
    if @format == :dat
      return to_dat
    end
    if @format == :csv
      return to_csv
    end
    if @format == :tex
      retval << generate_tex_table_header(@tex_page, @nosideways)
    end
    retval << format_table(@format)
    if @format == :tex
      retval << generate_tex_table_footer(@tex_page, @nosideways)
      retval.escape_for_latex!
    end
    return retval
  end

  private

  def determine_labels_col_width
    @labels_col_width = 0
    @table.each do |row|
      @labels_col_width = row.first.size if row.first.size > @labels_col_width
    end
    @labels_col_width += LabelsColSpacing
  end

  def determine_data_col_width
    @data_col_width = 0
    @table.each do |row|
      row.each_with_index do |cell, i|
        if i > 0
          @data_col_width = cell.size if cell.size > @data_col_width
        end
      end
    end
    @data_col_width += DataColSpacing
  end

  def to_csv
    rows = []
    @table.each do |values|
      values.each_with_index do |value, i|
        rows[i] ||= []
        rows[i] << value.to_s
      end
    end
    rows[0] = rows[0].map{|x| x.gsub(" ", "_")}
    retval = ""
    rows.each do |row|
      retval << row.join(",")
      retval << "\n"
    end
    return retval
  end

  def to_dat
    rows = []
    @table.each do |values|
      values.each_with_index do |value, i|
        rows[i] ||= []
        rows[i] << value.to_s
      end
    end
    rows[0] = rows[0].map{|x| x.gsub(" ", "_")}
    rows[0].unshift("#")
    retval = ""
    rows.each do |row|
      retval << row.join(" ")
      retval << "\n"
    end
    return retval
  end

  def format_table(output_format)
    retval = ""
    @table.each do |row|
      labels_col = true
      row.each do |col|
        (retval << sprintf(" & ")) if output_format == :tex && !labels_col
        width = labels_col ? @labels_col_width : @data_col_width
        retval << sprintf("%-#{width}s", col)
        labels_col = false
      end
      (retval << sprintf(" \\\\")) if output_format == :tex
      retval << sprintf("\n")
    end
    return retval
  end

  def generate_tex_table_header(tex_page, nosideways)
    retval = ""
    if (tex_page == :completedocument || tex_page == :first)
      retval << <<-endheader
\\documentclass[a4paper, 10pt]{article}
\\usepackage{booktabs}
\\usepackage{rotating}
\\pagestyle{empty}
\\begin{document}
      endheader
    elsif (tex_page == :middle || tex_page == :last)
      retval << "\\newpage\n"
    end
    if tex_page != :none
      retval << "\\begin{#{nosideways ? '' : 'sideways'}table}\n"
    end
    retval << "\\begin{tabular}{#{'l'*@nofcolumns}}\n"
    return retval
  end

  def generate_tex_table_footer(tex_page, nosideways)
    retval = ""
    retval << "\\end{tabular}\n"
    if tex_page != :none
      retval << "\\end{#{nosideways ? '' : 'sideways'}table}\n"
    end
    if (tex_page == :completedocument || tex_page == :last)
      retval << "\\end{document}\n"
    end
    return retval
  end

end

## DataCollection: data structure to contain the extracted data ##

class DataCollection

  OutOrder = [:title, :pipeline, :readset, :nofreads, :readlen, :tlen,
    :contreads, :contreads_percent, :minlen, :trovl, :introvl, :parts,
    :sga_index_d, :prefilter_time, :prefilter_space, :index_time,
    :index_space, :overlap_time, :overlap_space, :index_and_overlap_time,
    :index_and_overlap_space, :assembly_time, :assembly_space, :overall_time,
    :overall_space, :nofcontigs, :tcontigslen, :n50, :lngst]

  AutoCategory = {
    /time/ => :time,
    /space/ => :space,
    /mmap/ => :space,
    /alloc/ => :space,
    /vm_peak/ => :space,
    /hostname/ => :string,
    /compiler/ => :string,
    /indexname/ => :string,
    /lastcommit/ => :string,
    /percent/ => :percent,
    /hwm_gtsp_delta/ => :percent,
    /max\w+sort/ => :integer,
    /bits/ => :integer,
  }

  Category = {
    :title                      => :string,
    :minlen                     => :integer,
    :tlen                       => :millions,
    :readlen                    => :integer,
    :nofreads                   => :millions,
    :contreads                  => :thousands,
    :trovl                      => :millions,
    :introvl                    => :millions,
    :nofcontigs                 => :integer,
    :tcontigslen                => :longlength,
    :n50                        => :length,
    :lngst                      => :length,
    :parts                      => :integer,
    :sga_index_d                => :thousands,
    :prefixlength               => :integer,
    :sfxtotallength             => :millions,
    :sat                        => :string,
    :readset                    => :string,
  }

  Label = {
    :title                      => "",
    :minlen                     => "Min SPM length",
    :nofreads                   => "Number of reads",
    :tlen                       => "Total length",
    :readlen                    => "Read length",
    :contreads                  => "Contained reads",
    :contreads_percent          => "Contained %",
    :trovl                      => "Transitive SPM",
    :introvl                    => "Irreducible SPM",
    :nofcontigs                 => "Contigs",
    :tcontigslen                => "Tot contigs len",
    :n50                        => "Assembly N50",
    :lngst                      => "Longest contig",
    :parts                      => "Index parts",
    :sga_index_d                => "Sga index -d",
    :_time                      => "__phase_label__ time",
    :_t_percent                 => "__phase_label__ time (%)",
    :_gt_real_time              => " - gt:real",
    :_gt_cpu_time               => " - gt:CPU",
    :_gt_user_time              => " - gt:user",
    :_gt_sys_time               => " - gt:sys",
    :_real_time                 => " - real",
    :_cpu_time                  => " - CPU",
    :_user_time                 => " - user",
    :_sys_time                  => " - sys",
    :_real_cpu_time_diff        => " - real-CPU",
    :_gt_real_cpu_time_diff     => " - gt:real-CPU",
    :_space                     => "__phase_label__ space",
    :_s_percent                 => "__phase_label__ space (% peak)",
    :_gt_space                  => " - gt",
    :_hwm_space                 => " - vm_hwm",
    :_vm_peak                   => " - vm_peak",
    :_alloc                     => " - alloc",
    :_mmap                      => " - mmap",
    :_hwm_gtsp_delta            => " - hwm-gt (%)",
    :_indexname                 => "Indexname [__phase_short__]",
    :_hostname                  => "Hostname [__phase_short__]",
    :_compiler                  => "Compiler [__phase_short__]",
    :_lastcommit                => "Git commit [__phase_short__]",
    :_bits                      => "Ulong size (bits) [__phase_short__]",
    :_pipeline                  => "Pipeline [__phase_short__]",
    :indexname                  => "Indexname",
    :hostname                   => "Hostname",
    :compiler                   => "Compiler",
    :lastcommit                 => "Git commit",
    :pipeline                   => "Pipeline",
    :readset                    => "Readset",
    :bits                       => "Ulong size (bits)",
    :prefixlength               => "Sfx prefixlength",
    :maxinsertionsort           => "Sfx maxinssort",
    :maxbltriesort              => "Sfx maxbltsort",
    :maxcountingsort            => "Sfx maxcousort",
    :relevant                   => "Sfx relevant suf (%)",
    :suftabuint                 => "Sfx suftabuint",
    :memlimit                   => "Sfx memlimit",
    :sfxtotallength             => "Sfx totallength",
    :sat                        => "Sfx seq.ac.type",
 }

 UnitLabel = {
    :b                   => " (bytes)",
    :kb                  => " (kb)"  ,
    :Mb                  => " (Mb)"  ,
    :Gb                  => " (Gb)",
    :Tb                  => " (Tb)",
    :n                   => "",
    :K                   => " (K)",
    :M                   => " (M)",
    :bp                  => " (bp)",
    :kbp                 => " (kbp)",
    :Mbp                 => " (Mbp)",
    :Gbp                 => " (Gbp)",
    :s                   => " (s)",
    :hm                  => "",
  }

  RatioCategories = [:time, :space, :millions]

  NoValue = "-"

  def initialize
    @unit = {:millions => :M, :thousands => :K, :time => :s, :space => :Mb,
             :length => :kbp, :longlength => :Mbp}
    @extended = {}
    @hide = []
    @show = []
    @manual_cells = {}
    @space_variety = :hwm_or_gt
    @time_variety = nil
    @series       = []
    @outorder     = OutOrder.clone
    @manual_rows  = []
    @customlabels = {}
    @table        = []
    @keys         = []
    @index_and_overlap = false
    @show_keys = false
    @only_rows = []
    @ratios = []
    @ratiorows = []
    @main_phases = [:prefilter, :index, :overlap, :assembly, :overall]
    @phases = {
      :prefilter => {:order => 1, :label => "Prefilter",   :short => "P"},
      :index      => {:order => 2, :label => "Index",       :short => "I"},
      :overlap    => {:order => 3, :label => "Overlap",     :short => "O"},
      :assembly   => {:order => 4, :label => "Assembly",    :short => "A"},
      :overall    => {:order => 9, :label => "Overall",     :short => "all"},
    }
    @comb_phases = @main_phases.choose(2) + @main_phases.choose(3) +
      @main_phases.choose(4)
    @comb_phases.each do |phlist|
      key = phlist.map{|x|x.to_s}.join("_and_")
      order = @main_phases.index(phlist[0]) + 1
      label = phlist.map{|x|@phases[x][:label]}.join("/")
      short = phlist.map{|x|@phases[x][:short]}.join
      @phases[key] = {:order => order, :label => label, :short => short}
    end
    @phasewisekeys = [:indexname, :pipeline, :hostname, :compiler, :lastcommit,
      :bits]
  end

  attr_accessor :unit, :extended, :hide, :show, :manual_cells, :space_variety,
    :time_variety, :show_keys, :index_and_overlap, :only_rows, :manual_rows,
    :ratios, :ratiorows

  def create_series(name, title)
    title = name if @show_keys
    @series << [name, {:title => title}]
  end

  def add_value(name, key, value)
    if @series.assoc(name)[1][key].nil?
      @series.assoc(name)[1][key] = value
    end
  end

  def add_ratios
    pairs = []
    @ratiorows.map! do |str|
      str.to_sym
    end
    (0).step(@ratios.size-1,4) {|i| pairs << @ratios.slice(i,4)}
    pairs.each do |a, b, insert_after, ti|
      if !a.nil? and !b.nil?
        ratio_col = [:"ratio_#{a}_#{b}", {:ratio_column => true}]
        a_col = @series.assoc(a)
        b_col = @series.assoc(b)
        next if a_col.nil? || b_col.nil?
        if @ratiorows.empty?
          keys = a_col[1].keys & b_col[1].keys
          keys = keys.select{|k| RatioCategories.include?(category(k))}
        else
          keys = a_col[1].keys & b_col[1].keys & @ratiorows
        end
        keys.each do |key|
          a_value = a_col[1][key]
          next if !a_value.kind_of? Numeric
          b_value = b_col[1][key]
          next if !b_value.kind_of? Numeric || b_value == 0
          ratio_col[1][key] = a_value.to_f / b_value.to_f
        end
        if @show_keys
          ratio_col[1][:title] = ratio_col[0].to_s
        elsif ti == '-'
          autoti = "#{a_col[1][:title]}/#{b_col[1][:title]}"
          ratio_col[1][:title] = autoti if autoti != "/"
        else
          ratio_col[1][:title] = ti
        end
        insertion_point = @series.index(@series.assoc(insert_after)) + 1
        @series.insert(insertion_point, ratio_col)
      end
    end
  end

  def prepare!(rawmode)
    select_time_values(@time_variety)
    select_space_values(@space_variety)
    add_index_and_overlap_time
    add_index_and_overlap_space
    add_overall_time
    add_overall_space
    add_ratios if !@ratios.empty?
    if rawmode
      puts to_yaml
      exit
    end
    if !@only_rows.nil? && !@only_rows.empty?
      @outorder = @only_rows.map{|key|key.to_sym}
    end
    apply_manual_rows
    apply_manual_cells
    apply_extended
    @show.each {|key| append_row(key.to_sym) }
    apply_hide_option
    @outorder.uniq!
  end

  def tabulize
    prepare_first_column
    apply_space_time_extended
    maxwidth = 0
    @series.each do |n, values|
      @keys.each_with_index do |key, i|
        if values.has_key?(key)
          v = values[key]
          if values[:ratio_column]
            if key == :title || category(key) == :string
              str = format_value(v, :string)
            else
              str = format_value(v, :ratio)
            end
          elsif values[:"#{key}_string"]
            str = format_value(v, :string)
          else
            str = format_value(v, category(key))
          end
        else
          str = NoValue
        end
        str = str.to_s
        @table[i] << str
        maxwidth = str.size if str.size > maxwidth
      end
    end
    @colsize = maxwidth + 6
    return @table
  end

  private

  def category(key)
    return Category[key] if Category.has_key?(key)
    AutoCategory.each do |regex, categ|
      return categ if key.to_s =~ regex
    end
    return :string
  end

  def label(key)
    return @customlabels[key] if @customlabels.has_key?(key)
    return Label[key] if Label.has_key?(key)
    @phases.each_key do |phase|
      if key.to_s =~ /#{phase}(_.*)/
        subkey = $1.to_sym
        if Label.has_key?(subkey)
          template = Label[subkey].dup
          template.gsub!("__phase_label__", @phases[phase][:label])
          template.gsub!("__phase_short__", @phases[phase][:short])
          return template
        end
      end
    end
    return ""
  end

  def select_time_values(variety_to_use)
    if variety_to_use.to_s =~ /(.*)_or_gt_(.*)/
      select_time_values(:"gt_#{$1}")
      select_time_values(:"#{$1}")
    elsif variety_to_use.to_s =~ /gt_(.*)_or_(.*)/
      select_time_values(:"#{$1}")
      select_time_values(:"gt_#{$1}")
    else
      @series.each do |n, s|
        @main_phases.each do |ph|
          if s[:"#{ph}_#{variety_to_use}_time"]
            s[:"#{ph}_time"] = s[:"#{ph}_#{variety_to_use}_time"]
          end
        end
      end
    end
  end

  def select_space_values(variety_to_use)
    if variety_to_use == :hwm_or_gt
      select_space_values(:gt)
      select_space_values(:hwm)
    elsif variety_to_use == :gt_or_hwm
      select_space_values(:hwm)
      select_space_values(:gt)
    else
      @series.each do |n, s|
        @main_phases.each do |ph|
          if s[:"#{ph}_#{variety_to_use}_space"]
            s[:"#{ph}_space"] ||= s[:"#{ph}_#{variety_to_use}_space"]
          end
        end
      end
    end
  end

  def has?(key)
    @series.each do |n, values|
      if !values[key].nil?
        return true
      end
    end
    return false
  end

  def delete_row(key)
    pos = @keys.index(key)
    if pos
      @keys[pos] = nil
      @keys.compact!
      @table[pos] = nil
      @table.compact!
    end
  end

  def insert_row(key, after)
    if has?(key) && !@keys.include?(key)
      pos = @keys.index(after)
      if pos
        @keys.insert(pos + 1, key)
        @table.insert(pos + 1, [prepare_label(key)])
      end
    end
  end

  def append_row(key)
    if has?(key) && !@keys.include?(key)
      @keys << key
      @table << [prepare_label(key)]
    end
  end

  def apply_hide_option
    @hide.each do |key|
      if key =~ /(.*):(.*)/
        key = $1
        col = $2
        @series.assoc(col)[1][key.to_sym] = nil
      else
        @series.each do |n, values|
          values[key.to_sym] = nil
        end
      end
    end
  end

  def apply_manual_cells
    @manual_cells.each do |string|
      if string =~ /(.*):(.*):(.*)/
        key = $1
        col = $2
        value = $3
        @series.assoc(col)[1][key.to_sym] = value
        @series.assoc(col)[1][(key + "_string").to_sym] = true
      else
        STDERR.puts "Error: -setvalue #{string} not valid"
      end
    end
  end

  def apply_space_time_extended
    if @extended[:space]
      spkeys = %w{hwm_gtsp_delta vm_peak hwm_space gt_space mmap alloc}
      @main_phases.each do |ph|
        spkeys.each {|v| insert_row(:"#{ph}_#{v}", :"#{ph}_space")}
      end
      spkeys.reverse.each {|v| append_row(:"overlap_and_assembly_#{v}")}
    end
    if @extended[:s_percent] or @extended[:space]
      @series.each do |n, values|
        if !values[:ratio_column]
          time = values[:overall_space]
          if !time.nil?
            (@main_phases - [:overall]).each do |ph|
              v = values[:"#{ph}_space"]
              if !v.nil?
                values[:"#{ph}_s_percent"] = v.to_f * 100 / time if v
              end
            end
          end
        end
      end
      @main_phases.each {|ph| insert_row(:"#{ph}_s_percent", :"#{ph}_space")}
    end
    if @extended[:time]
      tk = %w{sys_time user_time real_cpu_time_diff cpu_time real_time
        gt_sys_time gt_user_time gt_real_cpu_time_diff gt_cpu_time
        gt_real_time}
      @main_phases.each do |ph|
        tk.each {|v| insert_row(:"#{ph}_#{v}", :"#{ph}_time")}
      end
      tk.reverse.each {|v| append_row(:"overlap_and_assembly_#{v}") }
    end
    if @extended[:t_percent] or @extended[:time]
      @series.each do |n, values|
        if !values[:ratio_column]
          time = values[:overall_time]
          if !time.nil?
            (@main_phases - [:overall]).each do |ph|
              v = values[:"#{ph}_time"]
              if !v.nil?
                values[:"#{ph}_t_percent"] = v.to_f * 100 / time if v
              end
            end
          end
        end
      end
      @main_phases.each {|ph| insert_row(:"#{ph}_t_percent", :"#{ph}_time")}
    end
  end

  def apply_extended
    extra_keys = @phasewisekeys.clone + [:prefixlength, :maxinsertionsort,
      :maxbltriesort, :maxcountingsort, :relevant, :suftabuint, :memlimit,
      :sat, :sfxtotallength]
    @phasewisekeys.each do |k|
      summarized = true
      @series.each do |n, values|
        vs = @phases.keys.map {|ph| values[:"#{ph}_#{k}"]}.compact!
        if vs.uniq.size == 1
          values[k] = vs[0]
        else
          summarized = false
        end
      end
      if !summarized
        extra_keys += @phases.keys.map {|ph| :"#{ph}_#{k}"}
      end
    end
    if @extended[:extra]
      extra_keys.each {|k| append_row(k)}
    end
  end

  def apply_manual_rows
    @manual_rows.each do |row|
      key = row.shift
      label = row.shift
      [row.size, @series.size].min.times do |i|
        @series[i][1][key] = row[i]
      end
      if label != "-" || label(key) == ""
        @customlabels[key] = label
      end
      if !@outorder.include?(key)
        @outorder << key
      end
    end
  end

  def prepare_label(key)
    if (@show_keys)
      if key == :title
        str = ""
      else
        str = key.to_s
      end
    else
      str = label(key)
      unit = UnitLabel[@unit[category(key)]]
      (str << unit) if !unit.nil?
    end
    return str
  end

  def prepare_first_column
    keys = []
    @series.each do |n, values|
      keys << values.keys
    end
    keys.uniq!
    keys.flatten!
    @outorder.each do |key|
      next if !has?(key)
      next if @keys.include?(key)
      if keys.include?(key)
        @keys << key
        @table << [prepare_label(key)]
      end
    end
    sizes = @table.map{|l| l[0].size}
    @firstcolsize = sizes.size > 0 ? sizes.max + 4 : 0
  end

  def format_number(number, divider, precision)
    return nil if !number
    number = number.to_f / divider
    number.with_precision(precision).to_s
  end

  def format_value(value, category)
    return NoValue if value.nil?
    case category
    when :time
      if @unit[:time] == :s
        return value.to_f.round.to_s
      elsif @unit[:time] == :hm
        return value.to_f.round.s2hm
      end
    when :space
      if @unit[:space] == :Mb
        return value.round.to_s
      elsif @unit[:space] == :Gb
        value = value.to_f.mem_convert(:Mb, :Gb)
        return value.to_f.with_precision(1).to_s
      end
    when :length
      if @unit[:length] == :kbp
        return format_number(value, 1_000, 1)
      elsif @unit[:length] == :bp
        return value.to_f.round.to_s
      end
    when :longlength
      if @unit[:longlength] == :Mbp
        return format_number(value, 1_000_000, 1)
      elsif @unit[:longlength] == :bp
        return value.to_f.round.to_s
      end
    when :millions
      if @unit[:millions] == :M
        return format_number(value, 1_000_000, 1)
      elsif @unit[:millions] == :n
        return value.to_f.round.to_s
      end
    when :thousands
      if @unit[:thousands] == :K
        return format_number(value, 1_000, 1)
      elsif @unit[:thousands] == :n
        return value.to_f.round.to_s
      end
    when :integer
      return value.to_f.round.to_s
    when :percent
      return value.to_f.with_precision(1).to_s
    when :string
      return value
    when :ratio
      precision = 2
      return sprintf("%.#{precision}f x", value.to_f.with_precision(precision))
    end
  end

  def overall_needed
    if @overall_needed.nil?
      @overall_needed = false
      @series.each do |n, values|
        skeys = 0
        tkeys = 0
        @main_phases.each do |ph|
          (skeys += 1) if values.has_key?(:"#{ph}_space")
          (tkeys += 1) if values.has_key?(:"#{ph}_time")
        end
        if skeys > 1 || tkeys > 1
          @overall_needed = true
          break
        end
      end
    end
    return @overall_needed
  end

  def add_index_and_overlap_space
    if @index_and_overlap
      @series.each do |n, values|
        spaces = []
        (spaces << values[:index_space]) if values[:index_space]
        (spaces << values[:overlap_space]) if values[:overlap_space]
        values[:index_and_overlap_space] = spaces.max
      end
    end
  end

  def add_index_and_overlap_time
    if @index_and_overlap
      @series.each do |n, values|
        time = 0
        (time += values[:index_time]) if values[:index_time]
        (time += values[:overlap_time]) if values[:overlap_time]
        time = nil if time == 0
        values[:index_and_overlap_time] = time
      end
    end
  end

  def add_overall_space
    if overall_needed
      @series.each do |n, values|
        spaces = []
        @main_phases.each do |ph|
          (spaces << values[:"#{ph}_space"]) if values[:"#{ph}_space"]
        end
        values[:overall_space] = spaces.max
      end
    end
  end

  def add_overall_time
    if overall_needed
      @series.each do |n, values|
        time = 0
        @main_phases.each do |ph|
          (time += values[:"#{ph}_time"]) if values[:"#{ph}_time"]
        end
        time = nil if time == 0
        values[:overall_time] = time
      end
    end
  end

end

class DataCollectionOptions
end

## ScriptRunner: parse arguments and coordinate execution ##

class ScriptRunner

  UsageMsg =<<-end_usage_msg
Collect information from the standard output of suffixerator/readjoiner and
some related programs in tabular form.

  Usage: #$0 [options] [<type>] <file> [[-<type>] <file>]+
end_usage_msg

  BasicUsageMsg=<<-end_usage_msg
#{UsageMsg}
  Use the -help option for more information.
end_usage_msg

  ExtendedUsageMsg=<<-end_usage_msg
#{UsageMsg}
  If two files for the same program are in phases order (e.g. I followed by O),
  they are grouped in a single column unless separated by --.

  An unique prefix of the filename is enough.

  Running time:
  -cpu    use cpu time = sum of user+sys [default]
  -real   use real time
  -user   use user time
  -tm     prefer time measurements by \"time\" over gt showtime [default]
  -gttm   prefer gt showtime over \"time\"
  -sec    express running time in seconds [default]
  -hm     express running time in hours and minutes
  -t%     additional show phase running time as percent of overall
  -t      output all available running times

  Space peak:
  -Mb     express space peaks in megabytes [default]
  -Gb     express space peaks in gigabytes
  -hwm    prefer vm_hwm over gt space peak [default]
  -gtsp   prefer gt space peak over vm_hwm
  -s%     additional show phase space peaks as percent of overall peak
  -s      output all available space peaks

  Further information:
  -m       use multipliers for lengths/counts (millions/thousands) [default]
  -n       display exact values of lenghts/counts
  -x       show extra information, such as hostname
  -io      show consolidate values for index + overlap

  Columns/Rows modifications:
  --       force start new column
  -ti xxx yyy ... --
           specify column titles
  -ratio key1 key2 key2 [title|-] --
           display ratio key2 / key1 after column key2
           either give ratio the default title (-) or the specified title
  -ratiorows key1 key2 ... --
           display ratio only for specified columns
  -hide key1 key2:column ... --
           hide specified rows (key1) or values (key2:column)
  -show key1 ... --
           additionally show specified rows (key1)
  -only key1 key2 .... --
           show only specified rows, in that order
  -setrow :symbol [label|-] value1 value2 ... --
           add specified row or manually set existing row to given values
  -setvalue key:column:value ... --
            manually set given values

  Tex table:
  -tex           output a latex table
  -no_page       output just the tabular content
  -first_page    skip page end
  -middle_page   skip page start and end
  -last_page     skip page end
  -no_sw         do not use sidewaystable

  Gnuplot data:
  -dat           output data for GnuPlot

  CSV:
  -csv           output data as comma delimited values

  File type override:
  -sfx       suffixerator
  -rdj[OA]   readjoiner
  -edn[OA]   edena
  -sga[IOA]  sga
  -sum       results file

  Further options:
  -keys    visualize row and column keys
  -raw     display the hash with the gathered data in yaml format and exit
  -clean   remove progress bars, \\r and trailing whitespaces from input files
           and exit (any other argument is considered to be a filename)
  -help    show this help message and exit

  Example usage:

  #$0 -tex -ti before_rebase after_rebase -- \\
       before_rebase.sfx.log before_rebase.rdj.log \\
       after_rebase.sfx.log after_rebase.rdj.log

  end_usage_msg

  LogTypes = ["-sfx", "-rdj", "-rdjO", "-rdjA", "-sum", "-ednO", "-ednA",
              "-sgaI", "-sgaO", "-sgaA", "-leap", "-pflt", "-uniq"]

  ParserClass = {
    "-sfx"  => OutputParser::Suffixerator,
    "-rdjO" => OutputParser::RdjOverlap,
    "-rdjA" => OutputParser::RdjAssembly,
    "-rdj"  => OutputParser::RdjBoth,
    "-sum"  => OutputParser::ResultsFile,
    "-ednO" => OutputParser::EdenaOverlap,
    "-ednA" => OutputParser::EdenaAssembly,
    "-sgaI" => OutputParser::SgaIndex,
    "-sgaO" => OutputParser::SgaOverlap,
    "-sgaA" => OutputParser::SgaAssembly,
    "-leap" => OutputParser::Leap,
    "-uniq" => OutputParser::Sequniq,
    "-pflt" => OutputParser::Prefilter,
  }

  SameColumn = {
    "-sfx"  => ["-rdj", "-rdjO"],
    "-rdjO" => ["-rdjA"],
    "-rdjA" => [],
    "-rdj"  => [],
    "-sum"  => [],
    "-ednO" => ["-ednA"],
    "-ednA" => [],
    "-sgaI" => ["-sgaO"],
    "-sgaO" => ["-sgaA"],
    "-sgaA" => [],
    "-leap" => [],
    "-uniq" => ["-sfx"],
    "-pflt" => ["-sfx", "-rdjO"],
  }

  Signatures =
  [
    # the first found in this order:
    ["prefixlength",                                "-sfx" ],
    ["gt readjoiner -overlap",                      "-rdjO"],
    ["gt readjoiner -assembly",                     "-rdjA"],
    ["gt readjoiner overlap",                       "-rdjO"],
    ["gt readjoiner assembly",                      "-rdjA"],
    ["gt readjoiner prefilter",                     "-pflt"],
    ["Computing overlaps",                          "-ednO"],
    ["loading transitively reduced overlap graph",  "-ednA"],
    ["Contigs elongations were stopped",            "-ednA"],
    ["calling mkqs",                                "-sgaI"],
    ["[sga::overlap]",                              "-sgaO"],
    ["sga assemble",                                "-sgaA"],
    ["Lngst",                                       "-sum" ],
    ["LEAP",                                        "-leap"],
    ["have been removed",                           "-uniq"],
    ["encseq2spm",                                  "-sfx" ],
  ]

  def parse!(args)
    parse_special_options!(args)
    parse_data_collection_options!(args)
    parse_output_options!(args)
    @titles = parse_stringarray!("-ti", args)
    @titles = nil if @titles.size == 0
    @data.hide += parse_stringarray!("-hide", args)
    @data.show = parse_stringarray!("-show", args)
    @data.only_rows = parse_stringarray!("-only", args)
    @data.ratios = parse_stringarray!("-ratio", args)
    @data.ratiorows = parse_stringarray!("-ratiorows", args)
    @data.manual_cells = parse_stringarray!("-setvalue", args)
    parse_setrow!(args)
    parse_series(args)
  end

  def initialize(arguments)
    @data = DataCollection.new
    @table = Table.new
    parse!(arguments)
    @series.each_with_index do |s, i|
      t = @titles ? @titles[i] : nil
      @data.create_series(i.to_s, t)
      s.each do |filename, logtype|
        log = IO.read(filename)
        parser = ParserClass[logtype].new(log)
        parser.results.each_pair do |key, value|
          @data.add_value(i.to_s, key, value)
        end
      end
    end
    @data.prepare!(@raw)
    @table.data = @data.tabulize
    puts @table.show
  end

  def new_series?(logtype, previous)
    return true if previous.nil? or logtype == previous
    return !SameColumn[previous].include?(logtype)
  end

  # delete progress bar lines and \r from transcripts
  def clean_script_transcript(filename)
    input = IO.read(filename)
    o = File.open(filename, "w")

    # remove \r:
    input.gsub!("\r\n", "\n")
    input.gsub!("\r", "\n")

    # delete progress bar and extra whitespace:
    input.each_line do |line|
      if line !~ /\d{1,3}\% \|[ \*]+\| .*ETA$/
        next if line =~ /cat: \/proc\/(\d+)\/status: No such file or directory/
        if line =~ / *\d{1,3}\% \|[ \*]+\| +(\d+:)?\d+:\d+ ?E?T?A?/
          line.gsub!(/ *\d{1,3}\% \|[ \*]+\| +(\d+:)?\d+:\d+ ?E?T?A?/, "")
        if line.gsub(/[\s:]/,"") != ""
          o.puts line.strip
        end
        else
          o.puts line.strip
        end
      end
    end

    o.close
  end

  def backup_transcript(filename, data)
    @backup ||= {}
    unless @backup[filename]
      b = File.open(filename + ".bak", "w")
      b.puts data
      b.close
      @backup[filename] = true
      STDERR.puts "The original transcript has been preserved"+
        "as #{filename}.bak"
    end
  end

  def clean_edena_transcript(filename)
    input = IO.read(filename)
    o = File.open(filename, "w")

    state = :normal
    previous_line = ""

    # delete progress bar and extra whitespace:
    input.lines.each_with_index do |line, lineno|
      case state
      when :normal
        if line =~ /^reading fasta entries/
          state = :reading_fasta_entries
        elsif line =~ /^Computing overlaps /
          state = :computing_overlaps
        elsif line =~ /^Removing transitive edges/
          state = :removing_transitive_edges
          next
        end
      when :reading_fasta_entries
        if line =~ /^fasta entries: /
          o.puts previous_line
          state = :normal
        elsif line =~ /^\d+/ or line.strip.empty?
          previous_line = line
          next
        else
          STDERR.puts "Warning: forced exit from "+
            "reading_fasta_entries state in #{filename}:#{lineno}"
          backup_transcript(filename, input)
          state = :normal
        end
      when :computing_overlaps
        if line =~ /^done/
          o.puts previous_line
          state = :normal
        elsif line =~ /^\d+/
          previous_line = line
          next
        else
          STDERR.puts "Warning: forced exit from "+
            "computing_overlaps state in #{filename}:#{lineno}"
          backup_transcript(filename, input)
          state = :normal
        end
      when :removing_transitive_edges
        if line =~ /^Removing transitive edges... done/
          state = :normal
        elsif line =~ /^Removing transitive edges...(\d+)/
          next
        else
          STDERR.puts "Warning: forced exit from "+
            "removing_transitive_edges state in #{filename}:#{lineno}"
          backup_transcript(filename, input)
          state = :normal
        end
      end
      o.puts line
    end

    o.close
  end

  def parse_special_options!(args)
    if args.size < 1
      puts BasicUsageMsg
      exit
    end
    if args.delete("-help")
      puts ExtendedUsageMsg
      exit
    end
    if args.delete("-clean")
      args.each do |filename|
        clean_script_transcript(filename)
        if ["-ednO"].include?(guess_logtype(IO.read(filename)))
          clean_edena_transcript(filename)
        end
      end
      exit
    end
  end

  def parse_data_collection_options!(args)
    @data.unit[:space] = :Mb if args.delete("-Mb")
    @data.unit[:space] = :Gb if args.delete("-Gb")
    @data.unit[:time] = :s if args.delete("-sec")
    @data.unit[:time] = :hm if args.delete("-hm")
    had_m = args.delete("-m")
    had_n = args.delete("-n")
    @data.unit[:length] = :kbp if had_m
    @data.unit[:length] = :bp if had_n
    @data.unit[:longlength] = :Mbp if had_m
    @data.unit[:longlength] = :bp if had_n
    @data.unit[:millions] = :M if had_m
    @data.unit[:millions] = :n if had_n
    @data.unit[:thousands] = :k if had_m
    @data.unit[:thousands] = :n if had_n
    v = "cpu" # default
    v = "cpu" if args.delete("-cpu")
    v = "real" if args.delete("-real")
    v = "user" if args.delete("-user")
    @data.time_variety = :"#{v}_or_gt_#{v}" # default
    @data.time_variety = :"#{v}_or_gt_#{v}" if args.delete("-tm")
    @data.time_variety = :"gt_#{v}_or_#{v}" if args.delete("-gttm")
    @data.space_variety = :hwm_or_gt if args.delete("-hwm")
    @data.space_variety = :gt_or_hwm if args.delete("-gtsp")
    @data.extended[:space] = (args.delete("-s"))
    @data.extended[:time] = (args.delete("-t"))
    @data.extended[:extra] = (args.delete("-x"))
    @data.extended[:t_percent] = (args.delete("-t%"))
    @data.extended[:s_percent] = (args.delete("-s%"))
    @data.show_keys = true if (args.delete("-keys"))
    @data.index_and_overlap = true if (args.delete("-io"))
    if @data.index_and_overlap
      @data.hide += [:index_time, :index_space, :overlap_time, :overlap_space]
    end
  end

  def parse_output_options!(args)
    @raw = (args.delete("-raw"))
    @table.format = :tex if args.delete("-tex")
    @table.format = :dat if args.delete("-dat")
    @table.format = :csv if args.delete("-csv")
    @table.tex_page = :first if (args.delete("-first_page"))
    @table.tex_page = :middle if (args.delete("-middle_page"))
    @table.tex_page = :last if (args.delete("-last_page"))
    @table.tex_page = :none if (args.delete("-no_page"))
    @table.nosideways = (args.delete("-nosw"))
  end

  def parse_stringarray!(starter, args)
    samode = false
    results = []
    args.each_with_index do |a, i|
      if !samode
        if a == starter
          samode = true
          args[i] = nil
        end
      else
        if a == "--"
          samode = false
        else
          results << a
        end
        args[i] = nil
      end
    end
    args.compact!
    return results
  end

  def parse_setrow!(args)
    all_setrow = parse_stringarray!("-setrow", args)
    @data.manual_rows = []
    all_setrow.each do |a|
      if a =~ /:(.*)/
        @data.manual_rows << [$1.to_sym]
      else
        if @data.manual_rows.empty?
          STDERR.puts "-setrow: missing \":<key>\""
          exit(-1)
        end
        @data.manual_rows.last << a
      end
    end
  end

  def guess_logtype(log)
    Signatures.each do |key, value|
      if log.index(key)
        logtype = value
        return logtype
      end
    end
    return nil
  end

  def parse_series(args)
    logtype = nil
    previous = nil
    @series = []
    args.each do |a|
      if LogTypes.include?(a)
        logtype = a
        next
      end
      if a == "--"
        @series << [] unless @series.size == 0 or @series.last.size == 0
        next
      end
      if !logtype
        if (!File.exists?(a))
          g = Dir.glob("#{a}*")
          if (g.size == 1)
            a = g[0]
          else
            STDERR.puts \
              "File \"#{a}\" not found (or "+
              "filename prefix not unique), skipped"
            next
          end
        end
        logtype = guess_logtype(IO.read(a))
      end
      if !logtype
        puts "Error: unknown logtype of #{a}, skipped"
        next
      end
      if new_series?(logtype, previous)
        @series << []
      end
      @series.last << [a, logtype]
      previous = logtype
      logtype = nil
    end
    if @series.empty?
      STDOUT.puts "Filenames list is empty, aborting..."
      exit(1)
    end
  end

end

if $0 == __FILE__
  ScriptRunner.new(ARGV)
end
