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
  io.stderr:write("Extract assemblegreedy program scaffold from GenomeTools " ..
                  "home directory gt_home.\n")
  os.exit(1)
end

if #arg == 1 then
  gt_home = arg[1]
else
  usage()
end

-- set up new project
name = "assemblegreedy"
p = Project:new(gt_home)
p:set_name(name)

-- add all dependencies
f = File:new("src/libgtcore/unused.h")
p:add(f)

f = File:new("src/libgtcore/undef.h")
p:add(f)

m = Module:new("src/libgtcore/xansi")
m:bare_includes()
p:add(m)

m = Module:new("src/libgtcore/countingsort")
m:bare_includes()
m:remove_include("error.h")
m:remove_function("get_int")
m:ma2xansi()
m:remove_unit_test()
p:add(m)

m = Module:new("src/libgtext/string_matching")
m:bare_includes()
m:remove_include("bittab.h")
m:remove_include("mathsupport.h")
m:remove_function("string_matching_bmh")
m:remove_function("string_matching_brute_force")
m:remove_function("string_matching_shift_and")
m:remove_function("store_first_match")
m:remove_function("store_match")
m:ma2xansi()
m:remove_unit_test()
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

-- add makefile
mf = Makefile:new(name, true)
mf:add_test("ruby -I. testsuite.rb")
p:set_makefile(mf)

-- add testsuite
p:add_stest()
testsuite = Testsuite:new("none_defined")
testsuite:add_file("example_fragments.fas", [[
>
aaaaacccccccccc
>
ggggggggggaaaaa
>
tttttttttaaaaac
>
ccccctttttggggg
]])
testsuite:add_file("example_fragments.assembly", [[
>Assembled sequence
tttttttttaaaaacccccccccctttttggggggggggaaaaa
]])
testsuite:add_test("assemblegreedy test", [[
  run "#{$bin}assemblegreedy #{$bin}example_fragments.fas"
  run "diff #{$last_stdout} #{$bin}example_fragments.assembly"
]])
p:add(testsuite)

-- add example program
prog = Program:new("assemblegreedy")
prog:add_include('<stdio.h>')
prog:add_include('<stdlib.h>')

prog:set_content([[
  if (argc != 2) {
    fprintf(stderr, "Usage: %s fragment_file\n", argv[0]);
    fprintf(stderr, "Assemble fragments from fragment_file in greedy "
                    "fashion\n");
    return EXIT_FAILURE;
  }

  /* XXX: your code gets called here */
]])

p:add(prog)

-- write tar
p:write_tar_file()
