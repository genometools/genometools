#!/usr/bin/env ruby
# what to do during travis tests

success = true
puts ENV['CC']
if ENV['CC'] == 'gcc' or ENV['GT_BITS'] == '32' or ENV['SYSTEM'] == 'Windows'
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
  if success and not ENV['SYSTEM'] == 'Windows' #win binaries won't run on linux
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

exit 0 if success
exit 1
