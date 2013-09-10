#!/usr/bin/env ruby
# what to do during travis tests

success = true
if ENV["CC"] == 'gcc'
  IO.popen([{"64bit"=>"yes"}, 'make', {:err=>[:child, :out]}]) do |io|
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
  IO.popen([{"64bit"=>"yes"}, 'make', 'test', {:err=>[:child, :out]}]) do |io|
    while (line = io.gets)
      print line
    end
  end
  success = $?.success?
end

if success
  exit 0
else
  exit 1
end
