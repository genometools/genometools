--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

if #arg == 2 then
  gt_home = arg[1]
  scorematrix_file = arg[2]
else
  io.stderr:write(string.format("Usage: %s gt_home scorematrix_file\n", arg[0]))
  io.stderr:write("Extract swalign program from GenomeTools home directory " ..
                  "gt_home.\n")
  os.exit(1)
end

-- returns C-function scorematrix_init() as string
local function generate_scorematrix_init()
  local sminit = {}
  local protein_alpha = gt.alpha_new_protein()
  local scorematrix = gt.scorematrix_read_protein(scorematrix_file)
  assert(protein_alpha:size() == scorematrix:get_dimension())
  sminit[#sminit + 1] = "static void scorematrix_init(int **scorematrix)"
  sminit[#sminit + 1] = "{"
  sminit[#sminit + 1] = "  assert(scorematrix);"
  for idx1 = 0, scorematrix:get_dimension()-1 do
    for idx2 = 0, scorematrix:get_dimension()-1 do
      sminit[#sminit + 1] = string.format("  scorematrix['%s']['%s'] = %d;",
                                          protein_alpha:decode(idx1),
                                          protein_alpha:decode(idx2),
                                          scorematrix:get_score(idx1, idx2))
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

-- add makefile
mf = Makefile:new(name)
p:set_makefile(mf)

-- add example program
prog = Program:new(name)
prog:add_include('<assert.h>')
prog:add_include('<limits.h>')
prog:add_include('<stdbool.h>')
prog:add_include('"array2dim.h"')
prog:add_include('"coordinate.h"')
prog:add_include('"minmax.h"')
prog:add_function(generate_scorematrix_init())
swtemp = File:new("src/libgtext/swalign.c")
prog:add_define("INDEL_SCORE", "-3")
prog:add_typedef(swtemp:get_typedef("DPentry"))
prog:add_function(swtemp:get_function("fillDPtable"))
prog:set_content([[
  int **scorematrix;
  DPentry **dptable;
  Coordinate alignment_end;
  const char *u, *v;
  unsigned long ulen, vlen;

  if (argc != 3) {
    fprintf(stderr, "Usage: %s protein_seq protein_seq\n", argv[0]);
    fprintf(stderr, "Align given protein sequences with Smith-Waterman.\n");
    return EXIT_FAILURE;
  }

  u = argv[1];
  v = argv[2];
  ulen = strlen(u);
  vlen = strlen(v);

  array2dim_calloc(scorematrix, CHAR_MAX, CHAR_MAX);
  array2dim_calloc(dptable, ulen+1, vlen+1);
  scorematrix_init(scorematrix);

  fillDPtable(dptable, u, ulen, v, vlen, (const int**) scorematrix,
              INDEL_SCORE, INDEL_SCORE, &alignment_end);

  /* XXX: include your code here... */

  array2dim_delete(dptable);
  array2dim_delete(scorematrix);

]])
p:add(prog)

-- write tar
p:write_tar_file()
