#!/usr/bin/env ruby
# add unique IDS to the fasta headers of file piped through stdin or file input
#

id = 1
ARGF.each_line do |line|
  if line.match(/^>\w+ (.*)$/)
    # print ">#{id.to_s(16)}#{id.to_s(16)}#{id.to_s(16)} #$1\n"
    # id += 1 + Math.log(id).to_i * Math.log(id).to_i * Math.log(id).to_i
    print ">0x#{id.to_s(16).rjust(4, '0')} #$1\n"
    id += 1
  else
    puts line
  end
end
