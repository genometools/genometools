--[[
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

-- testing the Lua bindings for the Encseq class

dnaseqfile = arg[1].."/lua_dnaseq.fas"
aaseqfile = arg[1].."/lua_aaseq.fas"
dseq1 = "agtccagctgtcagctagcgggcccgatgatatttt"
dseq2 = "gtgctgtac"
aaseq1 = "MVHFTAEEKAAVTSLWSKMNVEEAGGEALG"
aaseq2 = "KMNAVE"
idxsuffixes = {'esq','des','ssp','sds','al1'}

function run_test_num_seqs(es)
  assert(es:num_of_sequences() == 2)
end

function run_test_num_files(es)
  assert(es:num_of_files() == 1)
end

function run_test_descriptions(es)
  rval, err = pcall(GenomeTools_encseq.description, es, 2)
  assert(not rval)
  assert(string.find(err, "cannot exceed number of sequences"))
end

function run_test_total_length(es)
  assert(es:total_length() == 46)
end

function run_test_total_length_protein(es)
  assert(es:total_length() == 37)
end

function run_test_get_encoded_char(es, seq1, seq2)
  local a = es:alphabet()
  for i=1,seq1:len() do
    assert(a:decode(es:get_encoded_char(i-1, 0)) == seq1:sub(i,i))
    assert(es:get_decoded_char(i-1, 0) == seq1:sub(i,i))
  end
  for i=seq2:len(),1,-1 do
    assert(a:decode(es:get_encoded_char(seq2:len()-i, 1)) == seq2:sub(i,i))
    assert(es:get_decoded_char(seq2:len()-i, 1) == seq2:sub(i,i))
  end
  rval, err = pcall(GenomeTools_encseq.get_encoded_char, es, 100, 0)
  assert(not rval)
  assert(string.find(err, "cannot exceed"))
  rval, err = pcall(GenomeTools_encseq.get_encoded_char, es, 10, 6)
  assert(not rval)
  assert(string.find(err, "invalid readmode"))
end

function run_test_seq_startpos(es)
  assert(es:seqstartpos(0) == 0)
  assert(es:seqstartpos(1) == 37)
  rval, err = pcall(GenomeTools_encseq.seqstartpos, es, 2)
  assert(not rval)
  assert(string.find(err, "cannot exceed number of sequences"))
end

function run_test_seq_startpos_protein(es)
  assert(es:seqstartpos(0) == 0)
  assert(es:seqstartpos(1) == 31)
  rval, err = pcall(GenomeTools_encseq.seqstartpos, es, 2)
  assert(not rval)
  assert(string.find(err, "cannot exceed number of sequences"))
end

function run_test_seq_length(es)
  assert(es:seqlength(0) == 36)
  assert(es:seqlength(1) == 9)
  rval, err = pcall(GenomeTools_encseq.seqlength, es, 2)
  assert(not rval)
  assert(string.find(err, "cannot exceed number of sequences"))
end

function run_test_seq_length_protein(es)
  assert(es:seqlength(0) == 30)
  assert(es:seqlength(1) == 6)
  rval, err = pcall(GenomeTools_encseq.seqlength, es, 2)
  assert(not rval)
  assert(string.find(err, "cannot exceed number of sequences"))
end

function run_test_file_length(es)
  assert(es:effective_filelength(0) == 46)
  rval, err = pcall(GenomeTools_encseq.seqlength, es, 2)
  assert(not rval)
  assert(string.find(err, "cannot exceed"))
  rval, err = pcall(GenomeTools_encseq.get_encoded_char, es, 10, 6)
  assert(not rval)
  assert(string.find(err, "invalid readmode"))
end

function run_test_file_length_protein(es)
  assert(es:effective_filelength(0) == 37)
  rval, err = pcall(GenomeTools_encseq.seqlength, es, 2)
  assert(not rval)
  assert(string.find(err, "cannot exceed"))
end

function run_test_seq_substr_encoded(es, seq1, seq2)
  start = 3
  stop = 13
  res = es:extract_encoded(start, stop)
  a = es:alphabet()
  for i=start,stop do
    assert(a:decode(res[i-start+1]) == seq1:sub(i+1,i+1))
  end
  start = es:seqstartpos(1)
  stop = start + 4
  res = es:extract_encoded(start, stop)
  for i=start,stop do
    assert(a:decode(res[i-start+1]) == seq2:sub(i-start+1,i-start+1))
  end
  rval, err = pcall(GenomeTools_encseq.extract_encoded, es, 3, 1)
  assert(not rval)
  assert(string.find(err, "range endposition"))
  rval, err = pcall(GenomeTools_encseq.extract_encoded, es, 300, 500)
  assert(not rval)
  assert(string.find(err, "cannot exceed"))

end

function run_test_seq_substr_decoded(es, seq1, seq2)
  start = 3
  stop = 13
  res = es:extract_decoded(start, stop)
  a = es:alphabet()
  assert(res == seq1:sub(start+1,stop+1))
  start = es:seqstartpos(1)
  stop = start + 4
  res = es:extract_decoded(start, stop)
  assert(res == seq2:sub(1,5))
  rval, err = pcall(GenomeTools_encseq.extract_decoded, es, 3, 1)
  assert(not rval)
  assert(string.find(err, "range endposition"))
  rval, err = pcall(GenomeTools_encseq.extract_decoded, es, 300, 500)
  assert(not rval)
  assert(string.find(err, "cannot exceed"))
end

function run_test_seq_substr_sequential(es, seq1, seq2)
  start = 3
  stop = 13
  er = es:create_reader_with_readmode(0, start)
  a = es:alphabet()
  for i=start,stop do
    assert(a:decode(er:next_encoded_char()) == seq1:sub(i+1,i+1))
  end
  start = es:seqstartpos(1)
  stop = start + 4
  er = es:create_reader_with_readmode(0, start)
  for i=start,stop do
    encchar = a:decode(er:next_encoded_char())
    seqchar = seq2:sub(i-start+1, i-start+1)
    assert(encchar == seqchar)
  end
  rval, err = pcall(GenomeTools_encseq.create_reader_with_readmode, es, 0, 300)
  assert(not rval)
  assert(string.find(err, "cannot exceed"))
  rval, err = pcall(GenomeTools_encseq.create_reader_with_readmode, es, 7, 3)
  assert(not rval)
  assert(string.find(err, "invalid readmode"))
end

ee = gt.encseq_encoder_new()
ee:encode({dnaseqfile}, "dnaseqfile")
ee:encode({aaseqfile}, "aaseqfile")

el = gt.encseq_loader_new()
es = el:load("dnaseqfile")
run_test_descriptions(es)
run_test_num_seqs(es)
run_test_total_length(es)
run_test_num_files(es)
run_test_get_encoded_char(es, dseq1, dseq2)
run_test_seq_length(es)
run_test_seq_startpos(es)
run_test_file_length(es)
run_test_seq_substr_encoded(es, dseq1, dseq2)
run_test_seq_substr_decoded(es, dseq1, dseq2)
run_test_seq_substr_sequential(es, dseq1, dseq2)
a = gt.alphabet_new_dna()
eb = gt.encseq_builder_new(a)
eb:enable_multiseq_support()
eb:add_string(dseq1, "seq1")
eb:add_string(dseq2, "seq2")
es = eb:build()
run_test_descriptions(es)
run_test_num_seqs(es)
run_test_total_length(es)
run_test_num_files(es)
run_test_get_encoded_char(es, dseq1, dseq2)
run_test_seq_length(es)
run_test_seq_startpos(es)
run_test_file_length(es)
run_test_seq_substr_encoded(es, dseq1, dseq2)
run_test_seq_substr_decoded(es, dseq1, dseq2)
run_test_seq_substr_sequential(es, dseq1, dseq2)

es = el:load("aaseqfile")
run_test_descriptions(es)
run_test_num_seqs(es)
run_test_total_length_protein(es)
run_test_num_files(es)
run_test_get_encoded_char(es, aaseq1, aaseq2)
run_test_seq_length_protein(es)
run_test_seq_startpos_protein(es)
run_test_file_length_protein(es)
run_test_seq_substr_encoded(es, aaseq1, aaseq2)
run_test_seq_substr_decoded(es, aaseq1, aaseq2)
run_test_seq_substr_sequential(es, aaseq1, aaseq2)

a = gt.alphabet_new_protein()
eb = gt.encseq_builder_new(a)
eb:enable_multiseq_support()
eb:add_string(aaseq1, "seq1")
eb:add_string(aaseq2, "seq2")
es = eb:build()
run_test_descriptions(es)
run_test_num_seqs(es)
run_test_total_length_protein(es)
run_test_num_files(es)
run_test_get_encoded_char(es, aaseq1, aaseq2)
run_test_seq_length_protein(es)
run_test_seq_startpos_protein(es)
run_test_file_length_protein(es)
run_test_seq_substr_encoded(es, aaseq1, aaseq2)
run_test_seq_substr_decoded(es, aaseq1, aaseq2)
run_test_seq_substr_sequential(es, aaseq1, aaseq2)
