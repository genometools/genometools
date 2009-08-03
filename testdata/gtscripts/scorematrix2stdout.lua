--[[
  Copyright (c) 2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
  io.stderr:write("Parse score_matrix_file and show it on stdout.\n")
  os.exit(1)
end

if #arg == 1 then
  score_matrix_file = arg[1]
else
  usage()
end

blast_alphabet = gt.alphabet_new_empty()
blast_alphabet:add_mapping('A')
blast_alphabet:add_mapping('R')
blast_alphabet:add_mapping('N')
blast_alphabet:add_mapping('D')
blast_alphabet:add_mapping('C')
blast_alphabet:add_mapping('Q')
blast_alphabet:add_mapping('E')
blast_alphabet:add_mapping('G')
blast_alphabet:add_mapping('H')
blast_alphabet:add_mapping('I')
blast_alphabet:add_mapping('L')
blast_alphabet:add_mapping('K')
blast_alphabet:add_mapping('M')
blast_alphabet:add_mapping('F')
blast_alphabet:add_mapping('P')
blast_alphabet:add_mapping('S')
blast_alphabet:add_mapping('T')
blast_alphabet:add_mapping('W')
blast_alphabet:add_mapping('Y')
blast_alphabet:add_mapping('V')
blast_alphabet:add_mapping('B')
blast_alphabet:add_mapping('Z')
blast_alphabet:add_mapping('X')
blast_alphabet:add_wildcard('*')

blast_matrix = gt.score_matrix_new_read(score_matrix_file, blast_alphabet)
assert(blast_alphabet:size() == blast_matrix:get_dimension())

blast_matrix:show()
