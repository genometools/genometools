--[[
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
]]

require 'extractor' -- contains all the classes necessary for the extracting

function usage()
  io.stderr:write(string.format("Usage: %s gt_home\n", arg[0]))
  io.stderr:write("Extract linearalign program from GenomeTools home " ..
                  "directory gt_home.\n")
  os.exit(1)
end

if #arg == 1 then
  gt_home = arg[1]
else
  usage()
end

-- set up new project
name = "linearalign"
p = Project:new(gt_home)
p:set_name(name)

-- add all dependencies
f = File:new("src/libgtcore/unused.h")
p:add(f)

m = Module:new("src/libgtcore/xansi")
m:bare_includes()
p:add(m)

m = Module:new("src/libgtcore/dynalloc")
m:bare_includes()
m:ma2xansi()
p:add(m)

m = Module:new("src/libgtcore/array")
m:bare_includes()
m:remove_include("error.h")
m:remove_include("mathsupport.h")
m:remove_include("range.h")
m:remove_example()
m:remove_unit_test()
m:remove_function("array_iterate")
m:remove_function("iterate_test_func")
m:remove_function("iterate_fail_func")
m:ma2xansi()
p:add(m)

m = Module:new("src/libgtext/alignment")
m:bare_includes()
m:remove_include("error.h")
m:remove_unit_test()
m:ma2xansi()
p:add(m)

m = Module:new("src/libgtext/linearalign")
m:bare_includes()
m:ma2xansi()
p:add(m)

-- add makefile
mf = Makefile:new(name, true)
mf:add_test("ruby -I. testsuite.rb")
p:set_makefile(mf)

-- add testsuite
p:add_stest()
testsuite = Testsuite:new("none_defined")
testsuite:add_test("linearalign test 1", [[
  run "#{$bin}linearalign ab bbabcccbababcccb"
  run "diff #{$last_stdout} #{$bin}test1.out"
]])
testsuite:add_file("test1.out", [[
--ab------------
  ||            
bbabcccbababcccb
]])
testsuite:add_test("linearalign test 2", [[
  run "#{$bin}linearalign ccbaabac bbabcbabac"
  run "diff #{$last_stdout} #{$bin}test2.out"
]])
testsuite:add_file("test2.out", [[
ccba---abac
  ||   ||||
-bbabcbabac
]])
testsuite:add_test("linearalign test 3", [[
  run "#{$bin}linearalign bbababbbaba b"
  run "diff #{$last_stdout} #{$bin}test3.out"
]])
testsuite:add_file("test3.out", [[
bbababbbaba
|          
b----------
]])
p:add(testsuite)

-- add example program
prog = Program:new("main")
prog:add_include('<stdlib.h>')
prog:add_include('<string.h>')
prog:add_include('"linearalign.h"')

prog:set_content([[
  const char *u, *v;
  unsigned long ulen, vlen;
  Alignment *alignment;

  if (argc != 3) {
    fprintf(stderr, "Usage: %s seq_1 seq_2\n", argv[0]);
    fprintf(stderr, "Globally align seq_1 with seq_2 (with linear memory "
                    "consumption).\n");
    return EXIT_FAILURE;
  }

  u = argv[1];
  v = argv[2];
  ulen = strlen(u);
  vlen = strlen(v);

  alignment = linearalign(u, ulen, v, vlen);
  alignment_show(alignment, stdout);
  alignment_delete(alignment);

]])

p:add(prog)

-- write tar
p:write_tar_file()
