
# Patrick Maass, 2006

if $0 == __FILE__
  require 'stest'
  at_exit do
    OnError do exit 1 end
  end
end

def run_g(str, opts = {})
  cmd = {}
  prefix = ""
  if $arguments["valgrind"] == "yes"
    prefix = "valgrind --leak-check=yes "
  end
  opts[:stdout] = "std_out"
  opts[:stderr] = "std_err"
  run(prefix + str, opts)
  if $arguments["valgrind"] == "yes"
    run %{grep "ERROR SUMMARY: 0 errors" std_out std_err}
  end
end

def run_stest(fn, str, args = [], opts = {})
  tfn = "#{fn}.rb"
  File.open(tfn, "w") do |f|
    f.print(str)
  end
  run "ruby #{$runpath}/stest.rb #{tfn} #{args.join(" ")}", opts
end

Description "Hello, World"
Keywords "basics--apply--run"
Test do
  run_g "echo huhu"
  run "echo hello"
end

Name "Hello, World 2 and a long ending"
Keywords "basics--retval"
Test do
  run("grep foo barbaz", :retval => 2)
end

Descr "Unavailable executable"
Keywords "basics--no-exec"
Test do
  run_stest "tt", %{
    Name "No such exe"
    Keywords "test--no-exe"
    Test do
      run "no-such-exe foo"
    end
  }, [], :retval => 1
  grep($last_stdout, /problem: exec failed/)
end

Name "Unavailable executables"
Keywords "basics--multi-no-exec"
Test do
  run_stest "tt", %{
    Name "No such exe"
    Keywords "test--no-exe"
    Test do
      run "no-such-exe foo"
      run "no-such-exe-either"
    end
    Name "No such exe2"
    Keywords "test--no-exe2"
    Test do
      run "some-other bar"
    end
  }, [], :retval => 1
  grep($last_stdout, /problem: exec failed/)
  grep($last_stdout, /2:.*: error/)
end

Name "Arguments"
Keywords "basics--run-arguments"
Test do
  run "echo *"
  run "echo ______"
  run("echo", :arguments => [ "*" ])
end

Name "Environment"
Keywords "basics--environment"
Test do
  run("sh -c 'echo LL $TATA'", :env => { "TATA" => "foo" })
end

Name "Test Input Files"
Keywords "basics--write-file"
Test do
  rundir = $runpath
  File.open("tt.rb", "w") do |f|
    f.puts %{
      require 'stest'
      Name "TestTest"
      Test do
        run "echo hello"
      end
    }
  end
  run "ruby -I#{rundir} tt.rb"
end

Name "Test Error Log"
Keywords "basics--error-log"
Test do
  run_stest "tt", %{
    Name "Foo"
    Test do
      failtest "and this is why"
    end
  }, [], :retval => 1
  grep("stest_tt/test1/stest_error", /^and this is why:$/)
end

Name "Numbered output files"
Keywords "basics--output-files"
Test do
  run_stest "tt", %{
    Name "foo"
    Test do
      run "echo hello"
      run "echo huhu"
    end
  }
  run "test -r stest_tt/test1/stdout_1"
  run "test -r stest_tt/test1/stdout_2"
end

Name "Store commands"
Keywords "basics--store-commands"
Test do
  run_stest "tt", %{
    Name "foo"
    Test do
      run "echo hello"
    end
  }
  File.open("stest_tt/test1/run_1") do |f|
    l = f.readlines
    if l[0].chomp != "$CMD_PREFIX echo hello"
      failtest "invalid command stored"
    end
  end
end

Name "Test without keywords"
Keywords "basics--without-keywords"
Test do
  run_stest "tt", %{
    Name "foo"
    Test do
      run "echo"
    end
  }, [ "-keywords", "some" ]
  run "test ! -d stest_tt"
  run_stest "kk", %{
    Name "foo"
    Test do
      run "echo hello"
    end
  }, [ "-keywords", "not foo" ]
  run "test -r stest_kk/test1/run_1"
end

Name "Keyword Tokens"
Keywords "basics--keyword-parser"
Test do
  run_stest "tt", %{
    Name "foo"
    Keywords "x"
    Test do
      run "echo hello"
    end
    Name "bar"
    Keywords "notfoo foo"
    Test do
      run "echo huhu"
    end
  }, [ "-keywords", "notfoo" ]
  run "test ! -d stest_tt/test1"
  run "test -d stest_tt/test2"
end

Name "Test Selection"
Keywords "basics--selection"
Test do
  run_stest "tt", %{
    Name "s1"
    Test do
      run "echo foo"
    end
    Name "s2"
    Test do
      run "echo bar"
    end
    Name "s3"
    Test do
      run "echo baz"
    end
  }, %w{-select 2 3}
  grep("stest_tt/test2/stdout_1", "^bar$")
  grep("stest_tt/test3/stdout_1", "^baz$")
  run "test ! -d stest_tt/test1"
end

Name "Test Selection (name)"
Keywords "basics--name-selection"
Test do
  run_stest "tt", %{
    Name "foo"
    Test do
      run "echo hello"
    end
    Name "bar"
    Test do
      run "echo huhu"
    end
  }, [ "-name", "b.r" ]
  grep("stest_tt/test2/stdout_1", /^huhu$/)
  run "test ! -d stest_tt/test1"
end

Name "Inverse grep"
Keywords "utils--inverse-grep"
Test do
  run_stest "tt", %{
    Name "gt"
    Test do
      run "echo foo"
      grep($last_stdout, /foo/)
      grep($last_stdout, /bar/, true)
    end
  }
end

Name "Grep string"
Keywords "utils--string-grep"
Test do
  run_stest "tt", %{
    Name "sg"
    Test do
      run "echo foo"
      grep($last_stdout, ".")
    end
  }
end

Name "Multi-file grep"
Keywords "utils--multi-grep"
Test do
  run_stest "tt", %{
    Name "foo"
    Test do
      run "echo foo"
      grep([$last_stdout, $last_stderr], /^foo/)
      run "echo bar", :stdout => "s", :stderr => "e"
      grep(%w{e s}, /^bar/)
    end
  }
  run_stest "ft", %{
    Name "will-fail"
    Test do
      run "echo hello"
      grep([$last_stdout, $last_stderr], /nowhere/)
    end
  }, [], :retval => 1
  grep($last_stdout, /: failed/)
end

Name "Return Codes"
Keywords "basics--return-codes"
Test do
  run_stest "tt", %{
    Name "sg"
    Test do
      # will fail:
      run "echo foo", :retval => 1
    end
  }, [], :retval => 1
end

Name "Command Line Arguments"
Keywords "basics--command-line"
Test do
  run_stest "tt", %{
    Name "pa"
    Test do
      run "echo #\{$arguments["hello"]\}"
      run "echo #\{$arguments["flag"]\}"
    end
  }, %w{-hello hello world -flag}
  grep("stest_tt/test1/stdout_1", /^hello world$/)
  grep("stest_tt/test1/stdout_2", /^yes$/)
end

Name "Write-File"
Keywords "utils--write-file"
Test do
  run_stest "tt", %{
    Name "wf"
    Test do
      write_file "foo", %{
        from column 0
          indent
        end block
      }
    end
  }
  File.open("stest_tt/test1/foo", "r") do |f|
    l = f.readlines
    failtest "bad offset" if l[0] != "from column 0\n"
    failtest "line 1 != '  indent'" if l[1] != "  indent\n"
    failtest "line 2 != 'end block'" if l[2] != "end block\n"
    failtest "too many lines" if l.size != 3
  end
end

Name "Max Time"
Keywords "basics--max-time"
Test do
  run_stest "tt", %{
    Name "foo"
    Test do
      run "sleep 3", :maxtime => 1
    end
  }, [], :retval => 1
  grep $last_stdout, /did not finish on time/
end

Name "Range Selection"
Keywords "selection--ranges"
Test do
  run_stest "tt", %{
    Name "foo"
    Test {}
    Name "bar"
    Test {}
    Name "baz"
    Test {}
  }, %w{-select 2..3}
  grep $last_stdout, "bar"
  grep $last_stdout, "baz"
  grep $last_stdout, "foo", true
end

Name "Require/Provide"
Keywords "basics--require"
Test do
  run_stest "tt", %{
    Test.have_features do |str|
      if str =~ /^x/
        true
      else
        [ false, "skipped" ]
      end
    end
    Name "foo"
    Requires "x y v"
    Test do
      run "echo"
    end
    Name "bar"
    Requires "k y v"
    Test do
      run "echo"
    end
  }
  grep $last_stdout, /1:.*: ok/
end

