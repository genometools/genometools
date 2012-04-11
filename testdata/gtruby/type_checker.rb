#
# Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2012 Center for Bioinformatics, University of Hamburg
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

require 'gtruby'

if ARGV.size != 1 then
  STDERR.puts "Usage: #{$0} OBO_file"
  STDERR.puts "Test the type checker Ruby bindings using <OBO_file>."
  exit(1)
end

obofile = ARGV[0]

begin
  tc = GT::TypeChecker.new
rescue NotImplementedError => msg
  # expect exception
else
  raise TestFailedError
end

tc = GT::TypeCheckerBuiltin.new
raise unless tc.is_valid?("gene")
raise if tc.is_valid?("foo")

tc = GT::TypeCheckerOBO.new(obofile)
raise unless tc.is_valid?("gene")
raise unless tc.is_valid?("processed_transcript")
raise if tc.is_valid?("foo")

begin
  tc = GT::TypeCheckerOBO.new("fooo")
rescue GT::GTError => msg
  # expect exception
else
  raise TestFailedError
end
