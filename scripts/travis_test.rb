#!/usr/bin/env ruby
# what to do during travis tests

success = true
if ENV['CC'] == 'gcc' or ENV['GT_BITS'] == '32'
  env = {}
  if ENV['GT_BITS'] == '32'
    env['32bit' => 'yes']
  end
  IO.popen([ env,
            'make', '-j2',
            {:err=>[:child, :out]}]) do |io|
    while (line = io.gets)
      print line
    end
  end
  success = $?.success?
  if success
    IO.popen(["bin/gt", "-test", :err=>[:child, :out]]) do |io|
      while (line = io.gets)
        print line
      end
    end
    success = $?.success?
  end
else
  IO.popen(['make', '-j2', 'test',
            {:err=>[:child, :out]}]) do |io|
    while (line = io.gets)
      print line
      if m = line.match(/^\s*(\d+):\s.*: failed$/)
        errfiles = Dir.glob("testsuite/stest_testsuite/test#{m[1]}/stderr_*")
        errfiles.sort_by! do |path|
          path.match /\d+$/
          Integer($&)
        end
        puts "FAILING test, output of #{errfiles.last}"
        File.open(errfiles.last, 'r') do |file|
          file.each_line do |line|
            puts line
          end
        end
        puts "FAILING test, output of stest_error"
        File.open("testsuite/stest_testsuite/"+
                  "test#{m[1]}/stest_error", 'r') do |file|
          file.each_line do |line|
            puts line
          end
        end
      end
    end
  end
  success = $?.success?
end

if success
  exit 0
else
  exit 1
end
