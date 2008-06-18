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
  io.stderr:write("Extract consensus_sa program scaffold from GenomeTools " ..
                  "home directory gt_home.\n")
  os.exit(1)
end

if #arg == 1 then
  gt_home = arg[1]
else
  usage()
end

-- set up new project
name = "consensus_sa"
p = Project:new(gt_home)
p:set_name(name)

-- add all dependencies
f = File:new("src/libgtcore/fptr.h")
f:bare_includes()
p:add(f)

f = File:new("src/libgtcore/undef.h")
p:add(f)

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

m = Module:new("src/libgtcore/bittab")
m:bare_includes()
m:remove_include("mathsupport.h")
m:remove_unit_test()
m:remove_function("bittab_example")
m:ma2xansi()
m:fa2xansi()
p:add(m)

m = Module:new("src/libgtcore/error")
m:bare_includes()
m:remove_include("cstr.h")
m:ma2xansi()
p:add(m)

m = Module:new("src/libgtcore/range")
m:bare_includes()
m:remove_include("minmax.h")
m:remove_include("msort.h")
m:remove_include("safearith.h")
m:remove_function("range_compare_with_delta")
m:remove_function("generic_ranges_uniq")
m:remove_function("generic_ranges_uniq_in_place")
m:remove_function("ranges_uniq")
m:remove_function("ranges_uniq_count")
m:remove_function("ranges_uniq_in_place")
m:remove_function("ranges_uniq_in_place_count")
m:remove_function("range_offset")
m:remove_function("ranges_sort_by_length_stable")
m:remove_unit_test()
p:add(m)

m = Module:new("src/libgtexercise/sspliced_alignment")
m:bare_includes()
m:remove_include("cstr.h")
m:ma2xansi()
p:add(m)

m = Module:new("src/libgtexercise/sspliced_alignment_parsing")
m:bare_includes()
m:fa2xansi()
p:add(m)

m = Module:new("src/libgtcore/str")
m:bare_includes()
m:ma2xansi()
m:remove_include("cstr.h")
m:remove_include("ensure.h")
m:remove_include("error.h")
m:remove_include("genfile.h")
m:remove_function("str_clone")
m:remove_function("str_read_next_line_generic")
m:remove_function("str_unit_test")
p:add(m)

-- add makefile
mf = Makefile:new(name, true)
mf:add_test("ruby -I. testsuite.rb")
p:set_makefile(mf)

-- add testsuite
p:add_stest()
testsuite = Testsuite:new("none_defined")
testsuite:add_directory("testdata/consensus_sa", "testdata")
testsuite:add_content([[
i = 1
for infile in `ls #{$bin}testdata/*.in` do
  Name "consensus_sa test #{i}"
  Keywords "consensus_sa"
  Test do
    run "#{$bin}consensus_sa #{infile}"
    outfile = infile.gsub(/\.in$/, ".out")
    run "diff #{$last_stdout} #{outfile}"
  end
  i += 1
end
]])
--testsuite:
p:add(testsuite)

-- add example program
prog = Program:new("consensus_sa")
prog:add_include('<stdio.h>')
prog:add_include('<stdlib.h>')
prog:add_include('"array.h"')
prog:add_include('"error.h"')
prog:add_include('"fptr.h"')
prog:add_include('"sspliced_alignment.h"')
prog:add_include('"sspliced_alignment_parsing.h"')

prog:set_content([[
  Array *spliced_alignments;
  SSplicedAlignment *sa;
  unsigned long i;
  Error *err;
  int had_err = 0;

  if (argc != 2) {
    fprintf(stderr, "Usage: %s spliced_alignment_file\n", argv[0]);
    fprintf(stderr, "Read file containing spliced alingments, compute "
            "consensus spliced alignments,\nand print them to stdout.");
    return EXIT_FAILURE;
  }

  /* parse input file and store resuilts in the spliced alignment array */
  err = error_new();
  spliced_alignments = array_new(sizeof (SSplicedAlignment*));
  had_err = sspliced_alignment_parse(spliced_alignments, argv[1], err);

  if (!had_err) {
    /* sort spliced alignments */
    qsort(array_get_space(spliced_alignments), array_size(spliced_alignments),
          array_elem_size(spliced_alignments),
          (Compare) sspliced_alignment_compare_ptr);

    /* compute the consensus spliced alignments */
    /* XXX: your code gets called here */
  }

  /* free */
  for (i = 0; i < array_size(spliced_alignments); i++) {
    sa = *(SSplicedAlignment**) array_get(spliced_alignments, i);
    sspliced_alignment_delete(sa);
  }
  array_delete(spliced_alignments);

  /* report error, if necessary */
  if (had_err) {
    fprintf(stderr, "%s: error: %s\n", argv[0], error_get(err));
    error_delete(err);
    return EXIT_FAILURE;
  }
]])

p:add(prog)

-- write tar
p:write_tar_file()
