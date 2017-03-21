#!/usr/bin/env ruby

DISPLAY_OPTIONS = 
  "specify what information about the matches to display\n" +
  "alignment:    display alignment (possibly followed by =<number>\n" +
  "              to specify width of alignment columns)\n" +
  "cigar:        show cigar string representing alignment\n" +
  "polinfo:      display polishing information for displayed\n" +
  "              alignment\n" +
  "fstperquery:  output only the first found match per query\n" +
  "seed:         display the seed of the match, i.e. the length and\n" +
  "              the start position of the seed in both instances\n" +
  "failed_seed:  display the seed of the match that was extended,\n" +
  "              but failed (after extension) the filter conditions\n" +
  "seed_in_algn: display the seed in alignment\n" +
  "seqlength:    display length of sequences in which\n" +
  "              the two match-instances occur\n" +
  "evalue:       display evalue\n" +
  "s.seqdesc:    display sequence description of subject sequence\n" +
  "q.seqdesc:    display sequence description of query sequence\n" +
  "bitscore:     display bit score\n"

def keywords()
  extra = ["alignment","polinfo","fstperquery","failed_seed","seed_in_align"]
  kws = Array.new()
  idx = 0
  DISPLAY_OPTIONS.scan(/\n([a-z_\.]+):/) do |m|
    extra_flag = if extra.member?(m[0]) then "True" else "False" end
    kws.push([m[0],idx,extra_flag])
    idx += 1
  end
  return kws
end

puts keywords().sort {|a,b| a[0] <=> b[0]}.
                map{|s,idx,f| "{\"#{s}\", #{idx}, #{f}}"}.join(",\n")
