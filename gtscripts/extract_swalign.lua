--[[
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
  io.stderr:write(string.format("Usage: %s [-sol] gt_home score_matrix_file\n",
                                arg[0]))
  io.stderr:write("Extract swalign program from GenomeTools home directory " ..
                  "gt_home.\n")
  os.exit(1)
end

if #arg >= 2 then
  if string.match(arg[1], "^-sol") then
    solution = true
    table.remove(arg, 1)
  else
    solution = false
  end
  if #arg == 2 then
    gt_home = arg[1]
    score_matrix_file = arg[2]
  else
    usage()
  end
else
  usage()
end

-- returns C-function score_matrix_init() as string
local function generate_score_matrix_init()
  local sminit = {}
  local protein_alpha = gt.alpha_new_protein()
  local score_matrix = gt.score_matrix_new_read_protein(score_matrix_file)
  assert(protein_alpha:size() == score_matrix:get_dimension())
  sminit[#sminit + 1] = "static void score_matrix_init(int **score_matrix)"
  sminit[#sminit + 1] = "{"
  sminit[#sminit + 1] = "  assert(score_matrix);"
  for idx1 = 0, score_matrix:get_dimension()-1 do
    for idx2 = 0, score_matrix:get_dimension()-1 do
      sminit[#sminit + 1] = string.format("  score_matrix['%s']['%s'] = %d;",
                                          protein_alpha:decode(idx1),
                                          protein_alpha:decode(idx2),
                                          score_matrix:get_score(idx1, idx2))
    end
  end
  sminit[#sminit + 1] = "}"
  sminit[#sminit + 1] = ""
  return table.concat(sminit, "\n")
end

-- set up new project
name = "swalign"
p = Project:new(gt_home)
p:set_name(name)

-- add all dependencies
f = File:new("src/libgtcore/array2dim.h")
f:bare_includes()
f:remove_include("error.h")
f:remove_example()
f:ma2xansi()
p:add(f)

f = File:new("src/libgtcore/minmax.h")
p:add(f)

f = File:new("src/libgtcore/undef.h")
p:add(f)

f = File:new("src/libgtcore/unused.h")
p:add(f)

f = File:new("src/libgtext/coordinate.h")
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

-- add test files
f = File:new("testdata/gt_swalign_simple.out")
p:add(f)

-- add makefile
mf = Makefile:new(name)
mf:add_test([[
  ./swalign QSEIQANT SEQAN | diff - gt_swalign_simple.out
]])
p:set_makefile(mf)

-- add example program
prog = Program:new(name)
prog:add_include('<assert.h>')
prog:add_include('<limits.h>')
prog:add_include('<stdbool.h>')
prog:add_include('"alignment.h"')
prog:add_include('"array2dim.h"')
prog:add_include('"coordinate.h"')
prog:add_include('"minmax.h"')
prog:add_function(generate_score_matrix_init())
swtemp = File:new("src/libgtext/swalign.c")
prog:add_define("INDEL_SCORE", "-3")
prog:add_typedef(swtemp:get_typedef("DPentry"))
prog:add_function(swtemp:get_function("fillDPtable"))

if solution then
  prog:add_include('"undef.h"')
  prog:add_function(swtemp:get_function("traceback"))
  prog:add_function(swtemp:get_function("smith_waterman_align"))
end

swalign_start = [[
  int **score_matrix;
  DPentry **dptable;
  const char *u, *v;
  unsigned long ulen, vlen;
  Alignment *alignment = NULL;

  if (argc != 3) {
    fprintf(stderr, "Usage: %s protein_seq protein_seq\n", argv[0]);
    fprintf(stderr, "Align given protein sequences with Smith-Waterman.\n");
    return EXIT_FAILURE;
  }

  u = argv[1];
  v = argv[2];
  ulen = strlen(u);
  vlen = strlen(v);

  array2dim_calloc(score_matrix, CHAR_MAX, CHAR_MAX);
  array2dim_calloc(dptable, ulen+1, vlen+1);
  score_matrix_init(score_matrix);
  alignment = alignment_new();

]]

swalign_end = [[

  alignment_delete(alignment);
  array2dim_delete(dptable);
  array2dim_delete(score_matrix);

]]

if solution then
prog:set_content(swalign_start .. [[
  alignment = smith_waterman_align(u, v, u, v, ulen, vlen,
                                   (const int**) score_matrix,
                                   INDEL_SCORE, INDEL_SCORE);
  alignment_show(alignment, stdout);
  printf("\n");
]] .. swalign_end)
else
prog:set_content(swalign_start .. [[
  {
    Coordinate alignment_end;
    fillDPtable(dptable, u, ulen, v, vlen, (const int**) score_matrix,
                INDEL_SCORE, INDEL_SCORE, &alignment_end);

    /* XXX: include your code here... */
  }
]] .. swalign_end)
end

p:add(prog)

-- write tar
p:write_tar_file()
