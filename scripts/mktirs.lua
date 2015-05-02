#!/usr/bin/env gt
--[[
  Copyright (c) 2015 Sascha Steinbiss <sascha@steinbiss.name>

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

function usage()
  io.stderr:write("Creates artificial test sequences for TIRvish.\n")
  io.stderr:write(string.format("Usage: %s <prefix>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

prefix = arg[1]
totallen = 10000000
nof_tirs = 40
tirmin, tirmax = 200, 400
tirlmin, tirlmax = 1000, 6000
tsdmin, tsdmax = 6, 10
BASES = {'a','c','g','t'}
math.randomseed(os.time())

function revcomp(str)
  return string.gsub(str, "[acgt]",
                     {a = "t", c = "g", g = "c", t = "a"}):reverse()
end

function make_dna_of_len(l)
  out = {}
  for i = 1,l do
    out[i] = BASES[math.random(#BASES)]
  end
  return table.concat(out)
end

function make_tir()
  tirlen = tirmin + math.random(tirmax-tirmin)
  elem_len = tirlmin + math.random(tirlmax-tirlmin)
  tsdlen = tsdmin + math.random(tsdmax-tsdmin)
  tsd = make_dna_of_len(tsdlen)
  tir = make_dna_of_len(tirlen)
  elem_seq = tsd .. tir .. make_dna_of_len(elem_len-2*tirlen) ..
               revcomp(tir) .. tsd
  return elem_seq, string.len(elem_seq), tirlen, tir, revcomp(tir), tsdlen
end

base_seq = make_dna_of_len(totallen)
local annfile = assert(io.open(arg[1] .. ".gff3", "w+"))
seq = base_seq
for i = 1, nof_tirs do
  pos = math.random(string.len(seq)-1)
  elem_seq, elem_len, tirlen, tir, revtir, tsdlen = make_tir()
  annfile:write(prefix .."\t.\t" .. pos .. "\t" .. pos + elem_len -1 .. "\t"
                  .. ".\t+\t.\t")
  left = string.sub(seq, 1,pos)
  right = string.sub(seq, pos+1)
  seq = left .. elem_seq .. right
end
annfile:close()

local seqfile = assert(io.open(arg[1] .. ".fasta", "w+"))
seqfile:write(">" .. prefix .. "\n")
seqfile:write(seq)
seqfile:write("\n")
seqfile:close()
