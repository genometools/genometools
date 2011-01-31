#!/usr/bin/env ruby
class Sequence
	attr_reader :sequence, :length, :fasta_id
	attr_writer :sequence, :length, :fasta_id

	def initialize()
		@length = 0
		@sequence = ""
		@fasta_id = ""
	end

	def read_sequence(file)
		if !File.exists?(file)
			raise "The specified file #{file} does not exist."
		end
		File.open(file) do |line|
			#read first line
			str = line.gets.chomp
			if (str =~ /^>.*$/) then
				#FASTA file, read header properly
				@fasta_id = str
			else	
				#no FASTA file, but plain sequence
				@sequence = str
				@length = str.length
			end
			#read the rest of the file
			while nextline = line.gets
				nextline.chomp!
				@sequence = @sequence + nextline
				@length += nextline.length
			end
		end
	end
end

def revcomp(seq)
  seq.tr("aAcCgGtTnN","tTgGcCaAnN").reverse
end

class MSA
  attr_reader :seqs, :length, :cons

  def initialize(file)
    if !File.exists?(file)
      raise "The specified file #{file} does not exist."
    end
    @file = file
    @length = 0
    @seqs = []
    self.read
  end

  def read
    text = File.read(@file)
    text = text.split('>')
    text[1..text.length-1].each do |item|
      seq_components = item.split(/\n/)
      fasta_id=seq_components[0]
      sequence = seq_components[1..seq_components.length-1]
      seqObj = Sequence.new()
      seqObj.fasta_id = fasta_id
      seqObj.sequence = sequence.join
      seqObj.length = seqObj.sequence.length
      @seqs.push(seqObj)
    end
    @length = @seqs[0].length
  end
end

m = MSA.new(ARGV[0])
l = 0
m.seqs.each do |s|
  puts ">#{s.fasta_id}"
  puts "#{s.sequence}"
  l += s.length
end
m.seqs.reverse.each do |s|
  puts ">#{s.fasta_id}"
  puts "#{revcomp(s.sequence)}"
end

puts 
