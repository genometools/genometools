#
# Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2010 Center for Bioinformatics, University of Hamburg
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
require 'tempfile'

@dnaseqfile = Tempfile.new('dnaseq')
@aaseqfile = Tempfile.new('aaseq')
@dseq1 = "agtccagctgtcagctagcgggcccgatgatatttt"
@dseq2 = "gtgctgtac"
@dseq3 = "gtacagcac"
@dseq4 = "aaaatatcatcgggcccgctagctgacagctggact"
@dnaseqfile.write(">seq1\n#{@dseq1}\n")
@dnaseqfile.write(">seq2\n#{@dseq2}\n")
@aaseq1 = "MVHFTAEEKAAVTSLWSKMNVEEAGGEALG"
@aaseq2 = "KMNAVE"
@aaseqfile.write(">seq1\n#{@aaseq1}\n")
@aaseqfile.write(">seq2\n#{@aaseq2}\n")
@dnaseqfile.flush
@aaseqfile.flush
@idxsuffixes = ['esq','des','ssp','sds']

def create_es(indexname)
  GT::EncseqEncoder.new.encode([@dnaseqfile.path], indexname)
end

def create_es_protein(indexname)
  GT::EncseqEncoder.new.encode([@aaseqfile.path], indexname)
end
      
def create_mem
    a= GT::Alphabet.create_dna
    eb = GT::EncseqBuilder.create(a)
    eb.enable_description_support
    eb.enable_multiseq_support
    eb.add_string(@dseq1, 'seq1')
    eb.add_string(@dseq2, 'seq2')
    eb.build
end

def create_mem_protein
    a = GT::Alphabet.create_protein
    eb = GT::EncseqBuilder.create(a)
    eb.enable_description_support
    eb.enable_multiseq_support
    eb.add_string(@aaseq1, 'seq1')
    eb.add_string(@aaseq2, 'seq2')
    eb.build
end

def delete_idx(indexname)
    @idxsuffixes.each do |suf|
      File.unlink(indexname+"."+suf)
    end
end

val = create_es("foo")
@idxsuffixes.each do |suf|
  raise if !File.exists?("foo.#{suf}")
end
delete_idx("foo")

create_es("foo_mapped")
el = GT::EncseqLoader.new
es = el.load("foo_mapped")
raise if not es 
delete_idx("foo_mapped")

el = GT::EncseqLoader.new
begin
  el.load("foo_fail")
  raise
rescue
  # pass
end

def run_test_num_seqs(es)
  raise unless es.num_of_sequences == 2
  if es.alphabet.is_dna? then
    es.mirror
    raise unless es.num_of_sequences == 4
    es.unmirror
  end
end

def run_test_num_files(es)
  raise unless es.num_of_files == 1
  if es.alphabet.is_dna? then
    es.mirror
    raise unless es.num_of_files == 1
    es.unmirror
  end
end

def run_test_num_files_mem(es)
  raise unless es.num_of_files == 1
  if es.alphabet.is_dna? then
    es.mirror
    raise unless es.num_of_files == 1
    es.unmirror
  end
end

def run_test_descriptions(es)
  begin
    es.description(2)
    raise
  rescue
    #pass
  end
  raise unless es.description(0) == "seq1"
  raise unless es.description(1) == "seq2"
  if es.alphabet.is_dna? then
    es.mirror
    begin
      es.description(5)
      raise
    rescue
      #pass
    end
    raise unless es.description(3) == "seq1"
    raise unless es.description(2) == "seq2"
    es.unmirror
  end
end

def run_test_total_length(es)
  raise unless es.total_length == 46
  es.mirror
  raise unless es.total_length == 93
  es.unmirror
end

def run_test_total_length_protein(es)
  raise unless es.total_length == 37
end

def run_test_get_encoded_char(es, seq1, seq2, seq3=nil, seq4=nil)
  a = es.alphabet
  0.upto(seq1.length-1) do |i|
    encchar = es.get_encoded_char(i, GT::READMODE_FORWARD)
    raise unless a.decode(encchar) ==  seq1[i].chr
  end
  (seq2.length-1).downto(0) do |i|
    encchar = es.get_encoded_char(seq2.length-1-i, GT::READMODE_REVERSE)
    raise unless a.decode(encchar) == seq2[i].chr
  end
  if a.is_dna? then
    es.mirror
    0.upto(seq3.length-1) do |i|
      encchar = es.get_encoded_char(47+i, GT::READMODE_FORWARD)
      raise unless a.decode(encchar) ==  seq3[i].chr
    end
    (seq4.length-1).downto(0) do |i|
      encchar = es.get_encoded_char(seq4.length-1-i, GT::READMODE_REVERSE)
      raise unless a.decode(encchar) == seq4[i].chr
    end
    es.unmirror
  end
end

def run_test_seq_startpos(es)
  raise unless es.seqstartpos(0) == 0
  raise unless es.seqstartpos(1) == 37
  es.mirror
  raise unless es.seqstartpos(2) == 47
  raise unless es.seqstartpos(3) == 57
  es.unmirror
end

def run_test_seq_startpos_protein(es)
  raise unless es.seqstartpos(0) == 0
  raise unless es.seqstartpos(1) == 31
end

def run_test_seq_length(es)
  raise unless es.seqlength(0) == 36
  raise unless es.seqlength(1) == 9
  es.mirror
  raise unless es.seqlength(3) == 36
  raise unless es.seqlength(2) == 9
  es.unmirror
end

def run_test_file_length(es)
  raise "#{es.effective_filelength(0)} != 46" \
    unless es.effective_filelength(0) == 46
end

def run_test_seq_length_protein(es)
  raise unless es.seqlength(0) == 30
  raise unless es.seqlength(1) == 6
end

def run_test_file_length_protein(es)
  raise "#{es.effective_filelength(0)} != 37" \
    unless es.effective_filelength(0) == 37
end
  
def run_test_seq_substr_encoded(es, seq1, seq2, seq3=nil, seq4=nil)
  start = 3
  stop = 13
  res = es.extract_encoded(start, stop)
  a = es.alphabet()
  start.upto(stop) do |i|
    raise unless a.decode(res[i-start]) == seq1[i].chr
  end
  start = 0
  stop = 5
  spos = es.seqstartpos(1)
  res = es.extract_encoded(start+spos, stop+spos)
  start.upto(stop) do |i|
    raise unless a.decode(res[i-start]) == seq2[i].chr
  end
  if a.is_dna? then
    es.mirror
    spos = es.seqstartpos(2)
    start = 3
    stop = 8
    res = es.extract_encoded(start+spos, stop+spos)
    a = es.alphabet()
    start.upto(stop) do |i|
      raise unless a.decode(res[i-start]) == seq3[i].chr
    end
    start = 0
    stop = 5
    spos = es.seqstartpos(3)
    res = es.extract_encoded(start+spos, stop+spos)
    start.upto(stop) do |i|
      raise unless a.decode(res[i-start]) == seq4[i].chr
    end
    es.unmirror
  end
end

def run_test_seq_substr_plain(es, seq1, seq2, seq3=nil, seq4=nil)
  start = 3
  stop = 13
  raise unless es.extract_decoded(start, stop) == seq1[start..(stop)]
  start = 0
  stop = 5
  spos = es.seqstartpos(1)
  raise unless es.extract_decoded(spos+start, spos+stop) == seq2[start..(stop)]
  if es.alphabet.is_dna? then
    es.mirror
    spos = es.seqstartpos(2)
    start = 3
    stop = 8
    raise unless es.extract_decoded(start+spos, stop+spos) == seq3[start..(stop)]
    start = 0
    stop = 5
    spos = es.seqstartpos(3)
    raise unless es.extract_decoded(start+spos, stop+spos) == seq4[start..(stop)]
    es.unmirror
  end
end

def run_test_seq_substr_sequential(es, seq1, seq2, seq3=nil, seq4=nil)
  start = 3
  stop = 13
  er = es.create_reader_with_readmode(GT::READMODE_FORWARD, start)
  a = es.alphabet()
  start.upto(stop) do |i|
    raise unless a.decode(er.next_encoded_char) == seq1[i].chr
  end
  start = es.seqstartpos(1)
  stop = start + 5
  er = es.create_reader_with_readmode(GT::READMODE_FORWARD, start)
  0.upto(stop-start) do |i|
    encchar = a.decode(er.next_encoded_char)
    seqchar = seq2[i].chr
    raise unless encchar == seqchar
  end
  if a.is_dna? then
    es.mirror
    start = es.seqstartpos(2)
    stop = start + 5
    er = es.create_reader_with_readmode(GT::READMODE_FORWARD, start)
    0.upto(stop-start) do |i|
      encchar = a.decode(er.next_encoded_char)
      seqchar = seq3[i].chr
      raise unless encchar == seqchar
    end
    start = es.seqstartpos(3)
    stop = start + 5
    er = es.create_reader_with_readmode(GT::READMODE_FORWARD, start)
    0.upto(stop-start) do |i|
      encchar = a.decode(er.next_encoded_char)
      seqchar = seq4[i].chr
      raise unless encchar == seqchar
    end
    es.unmirror
  end
end

create_es("foo")
el = GT::EncseqLoader.new
es = el.load("foo")
run_test_descriptions(es)
run_test_get_encoded_char(es, @dseq1, @dseq2, @dseq3, @dseq4)
run_test_num_seqs(es)
run_test_num_files(es)
run_test_seq_length(es)
run_test_seq_startpos(es)
run_test_file_length(es)
run_test_seq_substr_encoded(es, @dseq1, @dseq2, @dseq3, @dseq4)
run_test_seq_substr_plain(es, @dseq1, @dseq2, @dseq3, @dseq4)
run_test_seq_substr_sequential(es, @dseq1, @dseq2, @dseq3, @dseq4)
run_test_total_length(es)
delete_idx("foo")
es = create_mem
run_test_descriptions(es)
run_test_get_encoded_char(es, @dseq1, @dseq2, @dseq3, @dseq4)
run_test_num_seqs(es)
run_test_num_files_mem(es)
run_test_seq_length(es)
run_test_seq_substr_encoded(es, @dseq1, @dseq2, @dseq3, @dseq4)
run_test_seq_substr_plain(es, @dseq1, @dseq2, @dseq3, @dseq4)
run_test_seq_substr_sequential(es, @dseq1, @dseq2, @dseq3, @dseq4)
run_test_total_length(es)

create_es_protein("foo")
el = GT::EncseqLoader.new
es = el.load("foo")
run_test_descriptions(es)
run_test_get_encoded_char(es,@aaseq1,@aaseq2)
run_test_num_seqs(es)
run_test_num_files(es)
run_test_seq_startpos_protein(es)
run_test_file_length_protein(es)
run_test_seq_length_protein(es)
run_test_seq_substr_encoded(es,@aaseq1,@aaseq2)
run_test_seq_substr_plain(es,@aaseq1,@aaseq2)
run_test_seq_substr_sequential(es,@aaseq1,@aaseq2)
run_test_total_length_protein(es)
delete_idx("foo")
es = create_mem_protein
run_test_descriptions(es)
run_test_get_encoded_char(es,@aaseq1,@aaseq2)
run_test_num_seqs(es)
run_test_num_files_mem(es)
run_test_seq_length_protein(es)
run_test_seq_substr_encoded(es,@aaseq1,@aaseq2)
run_test_seq_substr_plain(es,@aaseq1,@aaseq2)
run_test_seq_substr_sequential(es,@aaseq1,@aaseq2)
run_test_total_length_protein(es)
