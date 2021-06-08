#!/usr/bin/env ruby
# what to do during travis tests

success = true
opts = []
if !ENV["BUILDOPTS"].nil?
  opts = ENV["BUILDOPTS"].split(/\s+/)
end
call = ['make', 'test', 'testthreads=2', *opts, {:err=>[:child, :out]}]

IO.popen(call) do |io|
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

exit 0 if success
exit 1
