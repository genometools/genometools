#!/usr/bin/env ruby
# what to do during travis tests

success = true
if ENV["CC"] == 'gcc'
  print `make 64bit=yes 2>&1`
  success = $?.success?
  if success
    print `bin/gt -test 2>&1`
    success = $?.success?
  end
else
  print `make 64bit=yes test 2>&1`
  success = $?.success?
end

if success
  exit 0
else
  exit 1
end
