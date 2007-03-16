#!/usr/bin/env ruby

# Patrick Maass, 2006

require 'fileutils'

class TestFailed < Exception
end
class TestError < Exception
end

class Test
  @@id = 1
  @@prefix = "stest"
  @@failed = 0
  @@errors = 0
  @@feature_check = nil
  def self.prefix(str)
    @@prefix = str
  end
  def self.have_features(&bl)
    @@feature_check = bl
  end

  def initialize(str)
    @name = str
    @id = @@id
    @@id += 1
    @requirements = nil
  end

  def requires(str)
    @requirements = str
  end
  def keywords(str)
    @keywords ||= nil
    @did_run ||= nil
    raise "keyword reset" if @keywords
    raise "already run" if @did_run
    @keywords = {}
    str.split.each do |kw|
      if kw =~ /--/
        kws = kw.split(/--/)
        while kws.size > 0 do
          k = kws.join("--")
          @keywords[k] = true
          kws.pop
        end
      end
      @keywords[kw] = true
    end
  end

  def setup_dir
    begin
      ds = "#{@@prefix}/test#{@id}"
      if File.exists?(ds)
        FileUtils.rm_r(ds)
      end
      FileUtils.mkdir_p(ds)
      FileUtils.chdir(ds)
    rescue => e
      raise TestError, "could not create test dir '#{ds}'"
    end
  end

  def self.problems
    return @@failed, @@errors
  end


  def test(keywords, selection, name)
    raise "already run" if @did_run
    if name
      if not @name =~ name
        return
      end
    end
    if keywords
      @keywords ||= {}
      if not keywords.is_selected(@keywords)
        return
      end
    end
    if selection
      if not selection[@id]
        return
      end
    end
    if @requirements
      if @@feature_check
        ok, info = @@feature_check.call(@requirements)
        if not ok
          if info then
            STDOUT.printf "%3d: %-60s: #{info}\n", @id, @name
          end
          return
        end
      end
    end

    pid = fork do
      STDOUT.printf "%3d: %-60s: ", @id, @name
      STDOUT.flush
      e_info = nil
      e_code = 0
      begin
        setup_dir
        yield
      rescue TestFailed => e
        e_info = e
        e_code = 1
      rescue Interrupt => e
        exit 3
      rescue Exception => e
        e_info = e
        e_code = 2
      end
      if e_info
        i = "failed"
        if e_code != 1 then i = "error" end
        puts i
        puts "     [ problem: #{e_info.message}"
        puts "       in: #{Dir.pwd} ]"
        File.open("stest_error", "w") do |f|
          f.puts("Test #{@id} '#{@name}': #{i}:")
          f.puts("#{e_info.message}:")
          f.puts(e_info.backtrace)
        end
        exit e_code
      else
        puts "ok"
        exit 0
      end
    end
    begin
      p, s = Process::waitpid2(pid)
    rescue Interrupt => e
      puts "interrupted"
      exit 1
    end
    if s.exitstatus == 1
      @@failed += 1
    elsif s.exitstatus == 2
      @@errors += 1
    end
    @did_run = true
  end

end

def Test.exec(cmd, env, out_filename, err_filename, cmd_filename, max_t)
  pid = fork do
    env.each do |k,v|
      ENV[k] = v
    end
    begin
      File.open(cmd_filename, "w") do |f|
        env.each do |k, v|
          f.print "#{k}=#{v} "
        end
        # help for debugging failed tests:
        f.print "$CMD_PREFIX "
        f.puts(*cmd)
      end
      f_out = File.new(out_filename, "w")
      f_err = File.new(err_filename, "w")
      $stdout.reopen(f_out)
      $stderr.reopen(f_err)
    rescue Errno::ENOENT => e
      exit 90
    end
    begin
      Kernel.exec(*cmd)
    rescue Errno::ENOENT => e
      # XXX: better solution (signal/kill?)
      exit 91
    end
  end
  time_left = max_t
  sleep_step = 0.1
  s = false
  while time_left > 0
    sleep(sleep_step)
    time_left -= sleep_step
    p, s = Process.waitpid2(pid, Process::WNOHANG)
    break if s
  end
  if not s
    Process.kill("HUP", pid)
    sleep(sleep_step)
    p, s = Process.waitpid2(pid, Process::WNOHANG)
    raise TestError, "'#{cmd}' did not finish on time"
  end
  if not s.exited?
    if s.signaled?
      raise TestError, "'#{cmd}' received signal #{s.termsig}"
    else
      raise TestError, "'#{cmd}' stopped"
    end
  end
  if s.exitstatus == 90
    raise TestError, "cannot create output files for '#{cmd}'"
  end
  if s.exitstatus == 91
    raise TestError, "exec failed for '#{cmd}'"
  end
  if s.exitstatus == 92
    raise TestError, "exec error for '#{cmd}'"
  end
  return s.exitstatus
end

def Test.run(cmd, env, o_fn, r_fn, c_fn, rv, max_t)
  r = Test.exec(cmd, env, o_fn, r_fn, c_fn, max_t)
  if r != rv
    raise TestFailed, "unexpected return code: #{r} != #{rv}"
  end
end

class Keywords
  def initialize(kw_str)
    e, @root = OrNode.parse(kw_str, 0)
    if e != kw_str.size
      raise "extra tokens after '#{kw_str[0, e]}': '#{kw_str[e,kw_str.size]}'"
    end
  end
  def is_selected(valmap)
    return @root.is_selected(valmap)
  end

  class Node
    def self.skip_space(str, i)
      if i.class == NilClass
        raise "hu"
      end
      while str[i] and str[i] == ?\s
        i += 1
      end
      return i
    end
  end
  class VarNode < Node
    def self.parse(str, i)
      j = skip_space(str, i)
      if str[j] == ?(
        j += 1
        j, n = OrNode.parse(str, j)
        j = skip_space(str, j)
        if str[j] != ?)
          raise "missing ')'"
        end
        j += 1
        return j, n
      else
        s = j
        while str[j] and str[j] != ?\s and str[j] != ?( and str[j] != ?)
          j += 1
        end
        vn = str[s, j - s]
        if vn.size < 1 then raise "empty variable" end
        n = VarNode.new(vn)
        return j, n
      end
    end
    def is_selected(valset)
      if valset[@name] then return true end
      return false
    end
    def print
      puts "var: #{@name}"
    end
    private
    def initialize(name)
      @name = name
    end
  end
  class NotNode < Node
    def self.parse(str, i)
      is_neg = false
      j = skip_space(str, i)
      if str[j, 4] =~ /^not[ ()]?$/
        j += 3
        is_neg = true
      end
      j, n1 = VarNode.parse(str, j)
      if is_neg then n = NotNode.new(n1) else n = n1 end
      return j, n
    end
    def is_selected(valset)
      return !@n.is_selected(valset)
    end
    def print
      puts "not:"
      @n.print
    end
    private
    def initialize(n) @n = n end
  end
  class AndNode < Node
    def self.parse(str, i)
      j = skip_space(str, i)
      j, n1 = NotNode.parse(str, j)
      j = skip_space(str, j)
      if str[j, 4] =~ /^and[ ()]?$/
        j += 3
        j, n2 = AndNode.parse(str, j)
        n = self.new(n1, n2)
        return j, n
      end
      return j, n1
    end
    def is_selected(valset)
      return @n1.is_selected(valset) && @n2.is_selected(valset)
    end
    def print
      puts "A["
      @n1.print
      puts "and"
      @n2.print
      puts "]"
    end
    private
    def initialize(n1, n2) @n1 = n1; @n2 = n2 end
  end
  class OrNode < AndNode
    def self.parse(str, i)
      j = skip_space(str, i)
      j, n1 = AndNode.parse(str, j)
      j = skip_space(str, j)
      if str[j, 3] =~ /^or[ ()]?$/
        j += 2
        j, n2 = OrNode.parse(str, j)
        n = self.new(n1, n2)
        return j, n
      else
        return j, n1
      end
    end
    def is_selected(valset)
      return @n1.is_selected(valset) || @n2.is_selected(valset)
    end
    def print
      puts "O["
      @n1.print
      puts "or"
      @n2.print
      puts "]"
    end
  end

end

# Interface to Test class and helpers
# we use global variables for a simple user interface

$arguments = {}
$curr_test = nil
$run_no = 1

def run(str, opts = {})
  rv = 0
  if opts[:retval]
    rv = opts[:retval].to_i
  end
  o_fn = opts[:stdout]
  if not o_fn then o_fn = "stdout_#{$run_no}" end
  r_fn = opts[:stderr]
  if not r_fn then r_fn = "stderr_#{$run_no}" end
  c_fn = "run_#{$run_no}"
  $run_no += 1
  cmd = [ str ]
  if opts[:arguments]
    cmd += opts[:arguments]
  end
  env = opts[:env]
  env ||= {}
  max_t = opts[:maxtime]
  max_t ||= 10
  Test.run(cmd, env, o_fn, r_fn, c_fn, rv, max_t)
  $last_stdout = o_fn
  $last_stderr = r_fn
end

def failtest(msg)
  raise TestFailed, msg
end

def grep(f, patt, inverse=false)
  if not patt.respond_to?(:source)
    patt = Regexp.new(patt)
  end
  hit = nil
  if f.respond_to?(:to_str)
    files = [ f ]
  else
    files = f
  end
  files.each do |fn|
    File.open(fn) do |f|
      m = f.grep patt
      hit = m && m.size > 0
    end
    if hit then break end
  end
  p = "missing"
  if inverse then
    hit = !hit
    p = "unexpected"
  end
  if not hit
    raise TestFailed, "#{p} /#{patt.source}/ in #{files.join(" / ")}"
  end
end

def write_file(fn, str)
  File.open(fn, "w") do |f|
    # 1. remove first and last line
    # 2. remove offset
    a = str.split("\n")
    e = a.size - 1
    a[1] =~ /^( *)/
    off = $1.size
    a.each_with_index do |l, i|
      next if i == 0
      next if i == e
      f.puts "#{l[off, l.size]}"
    end
  end
end

def Name(str)
  t = Test.new(str)
  $curr_test = t
end

alias Description Name
alias Descr Name

def Keywords(str)
  t = $curr_test
  t.keywords(str)
end

def Requires(str)
  t = $curr_test
  t.requires(str)
end

def Test(&block)
  t = $curr_test
  t.test($arguments["keywords"], $arguments["select"],
         $arguments["name"], &block)
end

def OnError
  f, e = Test.problems
  if f > 0 || e > 0
    yield(f, e)
  end
end

# Set defaults

Test.prefix "stest_#{File.basename($0, ".rb")}"

$runpath = Dir.pwd

# Set arguments from command line

testfile = nil

if ARGV.size > 0
  i = 0
  while i < ARGV.size do
    a = ARGV[i]
    if a[0] == ?- then
      v = ARGV[i + 1]
      k = a[1..a.size]
      $arguments[k] ||= ""
      av = []
      while ARGV[i + 1] and ARGV[i + 1][0] != ?-
        av.push(ARGV[i + 1])
        i += 1
      end
      $arguments[k] = if av.size > 0 then av.join(" ") else "yes" end
    elsif $0 == __FILE__ and not testfile
      testfile = a
      Test.prefix "stest_#{File.basename(testfile, ".rb")}"
    else
      raise "invalid argument '#{a}'"
    end
    i += 1
  end
end

# Handle built-in options

if $arguments["keywords"]
  keywords = $arguments["keywords"]
  $arguments["keywords"] = Keywords.new(keywords)
end
if $arguments["select"]
  sel = $arguments["select"]
  $arguments["select"] = {}
  if sel =~ /(\d+)\.\.(\d+)/
    Range.new(Integer($1), Integer($2)).each do |id|
      $arguments["select"][Integer(id)] = true
    end
  else
    sel.split.each do |id|
      $arguments["select"][Integer(id)] = true
    end
  end
end
if $arguments["name"]
  $arguments["name"] = Regexp.new($arguments["name"])
end

if $0 == __FILE__
  if not testfile
    raise "missing test file argument"
  end
  require testfile
  OnError do exit(1) end
end

