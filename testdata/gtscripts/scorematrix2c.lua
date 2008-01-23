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

-- testing the Lua bindings for the ScoreMatrix class

function usage()
  io.stderr:write(string.format("Usage: %s score_matrix_file\n", arg[0]))
  io.stderr:write("Parse protein score_matrix_file and show it as C program.\n")
  os.exit(1)
end

if #arg == 1 then
  score_matrix_file = arg[1]
else
  usage()
end

protein_alpha = gt.alpha_new_protein()
score_matrix = gt.score_matrix_new_read_protein(score_matrix_file)
assert(protein_alpha:size() == score_matrix:get_dimension())

print("#include <limits.h>\n")
print("int main(int argc, char *argv[])")
print("{")
print("  int score_matrix[CHAR_MAX][CHAR_MAX] = {{0}};")
for idx1 = 0, score_matrix:get_dimension()-1 do
  for idx2 = 0, score_matrix:get_dimension()-1 do
    print(string.format("  score_matrix['%s']['%s'] = %d;",
          protein_alpha:decode(idx1), protein_alpha:decode(idx2),
          score_matrix:get_score(idx1, idx2)))
  end
end
print("  return 0;")
print("}")
