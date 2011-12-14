#!/usr/bin/env ruby
#
# Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
# Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg
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

helpmessage = <<-endhelpmessage
Align a set of contigs against a genomic sequence or simulation template
and output a list of misassembled contigs.

This script requires vmatch (http://www.vmatch.de/).

  Usage: #$0 <genomefilename> <contigsfilename> [maxedist_percent]

If an index <genomefilename>.vmatch already exists, it is used,
otherwise it is constructed using mkvtree.

By default only exact matches are considered correct. To accept approximate
matches use 0 < maxedist_percent < 100.
endhelpmessage

def get_template_index(templatefilename)
  indexname = "#{templatefilename}.vmatch"
  # the index is constructed only if necessary
  unless File.exists?("#{indexname}.suf")
    puts "# index #{indexname} not found, running mkvtree..."
    cmd = "mkvtree -dna -pl -allout -db #{templatefilename} "+
          "-indexname #{indexname}"
    puts cmd
    `#{cmd}`
  else
    puts "# using existing index #{indexname}"
  end
  return indexname
end

def run_vmatch(indexname, contigsfilename, maxedist_percent)
  if maxedist_percent
    edistopt = "-e #{maxedist_percent}b "
    desc = "approx"
  else
    edistopt = ""
    desc = "exact"
  end
  puts "# running vmatch (#{desc} matching)..."
  cmd = "vmatch -complete -d -p #{edistopt}-q #{contigsfilename} #{indexname}"
  puts cmd
  vmatch_output = `#{cmd}`
  # output to file
  filename = "#{indexname}.#{desc}.last_stdout"
  puts "# vmatch: done, transcript: #{filename}"
  outfile = File.open(filename, "w")
  outfile.puts vmatch_output
  outfile.close
  return vmatch_output
end

def parse_vmatch_output(vmatch_output)
  puts "# parsing vmatch output..."
  list = vmatch_output.scan(/[DP]\s+\d+\s+(\d+)/)
  list.flatten!
  list.map!{|x|Integer(x)}
  list.sort!
  list.uniq!
  return list
end

def get_nof_contigs(contigsfilename)
  cmd = "gt seqstat -contigs #{contigsfilename}"
  puts cmd
  stats = `#{cmd}`
  stats =~ /number of contigs:\s+(\d+)/
  nof_contigs = Integer($1)
  puts "# number of contigs: #{nof_contigs}"
  return nof_contigs
end

def find_not_in_list(list, nof_contigs, desc)
  puts "# create list of #{desc} contigs..."
  listcopy = list.dup
  not_in_list = []
  i = 0
  last = false
  while(!last) do
    in_list = listcopy.shift
    last = in_list.nil?
    in_list = nof_contigs if last
    i.upto(in_list - 1) do |x|
      not_in_list << x
    end
    i = in_list + 1
  end
  return not_in_list
end

def extract_from_fasta(infile, outfile, seqnums)
  seqnum = -1
  maxseqnum = seqnums.max
  infile.each do |line|
    seqnum += 1 if line[0] == ?>
    outfile.puts line if seqnums.include?(seqnum)
    break if seqnum > maxseqnum
  end
end

def output_contigs(contigsfilename, list, desc)
  filename = "#{contigsfilename}.#{desc}"
  puts "# writing to file: #{filename}"
  outfile = File.open(filename, "w")
  extract_from_fasta(IO.read(contigsfilename), outfile, list)
  outfile.close
  return filename
end

def correct_seqnums!(relative_seqnums, absolute_seqnums)
  relative_seqnums.map!{|i| absolute_seqnums[i]}
end

def run(templatefilename, contigsfilename, maxedist_percent)
  indexname = get_template_index(templatefilename)
  nof_contigs = get_nof_contigs(contigsfilename)
  vmatch_exact_output = run_vmatch(indexname, contigsfilename, nil)
  exact = parse_vmatch_output(vmatch_exact_output)
  nof_exact = exact.size
  puts "# number of exact contigs: #{nof_exact} (%.2f %%)" %
    ((nof_exact.to_f / nof_contigs) * 100)
  if nof_exact == nof_contigs
    puts "# all contigs are exact"
    exit
  end
  inexact = find_not_in_list(exact, nof_contigs, "inexact")
  nof_inexact = inexact.size
  puts "# number of inexact contigs: #{nof_inexact} (%.2f %%)" %
    ((nof_inexact.to_f / nof_contigs) * 100)
  inexact_contigsfilename =
    output_contigs(contigsfilename, inexact, "inexact")
  vmatch_approx_output = run_vmatch(indexname, inexact_contigsfilename,
                                    maxedist_percent)
  approx = parse_vmatch_output(vmatch_approx_output)
  nof_approx = approx.size
  puts "# number of approx contigs: #{nof_approx} (%.2f %%)" %
    ((nof_approx.to_f / nof_contigs) * 100)
  correct_seqnums!(approx, inexact)
  valid = approx + exact
  valid.sort!
  valid.uniq!
  nof_valid = valid.size
  puts "# number of valid contigs: #{nof_valid} (%.2f %%)" %
    ((nof_valid.to_f / nof_contigs) * 100)
  if nof_valid == nof_contigs
    puts "# all contigs are valid"
    exit
  end
  invalid = find_not_in_list(valid, nof_contigs, "invalid")
  nof_invalid = invalid.size
  puts "# number of invalid contigs: #{nof_invalid} (%.2f %%)" %
    ((nof_invalid.to_f / nof_contigs) * 100)
  invalid_contigsfilename =
    output_contigs(contigsfilename, invalid, "invalid")
end

# main:

if (ARGV.size < 2 || ARGV.size > 3)
  puts helpmessage
  exit
else
  run(ARGV[0], ARGV[1], ARGV[2])
end
