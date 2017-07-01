#!/usr/bin/env ruby

class String
  def dot2us
    return self.gsub(/\./,"_").gsub(/\s/,"")
  end
  def format_enum_value
    return "Gt_" + self.dot2us.capitalize + "_display"
  end
end

def keywords(display_options)
  kws = Array.new()
  idx = 0
  max_display_flag_length = 0
  display_options.each do |arg,helpline|
    incolumn = if helpline.match(/^display /) then "true" else "false" end
    kws.push([arg,idx,incolumn])
    idx += 1
    if max_display_flag_length < arg.length
      max_display_flag_length = arg.length
    end
  end
  return kws, max_display_flag_length
end

def indent(longest,arg)
  return " " * (longest - arg.length + 1)
end

def format(longest,helpline)
  len = longest + 3
  out = Array.new()
  helpline.split(/\s/).each do |w|
    if len + w.length <= 58
      out.push(w)
      len += w.length
    else
      out.push("\\n\"\n" + " " * 9 + "\"" + " " * (longest+2) + w)
      len = longest + 1 + w.length
    end
  end
  return out.join(" ").gsub(/ \\n/,"\\n")
end

TRACE_DELTA = 50
ALIGNMENT_WIDTH = 60

# The following array defines the keywords which can be used as arguments
# to option -outfmt, for each Keyword, a helpline is added.

display_options = [
  ["alignment",   "show alignment (possibly followed by =<number> to " +
                  "specify width of alignment columns, default is " +
                  "#{$ALIGNMENT_WDITH})"],
  ["seed_in_algn","mark the seed in alignment"],
  ["polinfo",     "add polishing information for shown alignment"],
  ["seed",        "abbreviation for seed.len seed.s seed.q"],
  ["failed_seed", "show the coordinates of a seed extension, which does not " +
                  "satisfy the filter conditions"],
  ["fstperquery", "output only the first found match per query"],
  ["tabsep",      "separate columns by tabulator, instead of space as default"],
  ["blast",       "output matches in blast format 7 (tabular with comment " +
                  "lines; instead of gap opens, indels are displayed)"],
  ["gfa2",        "output matches in gfa2 format"],
  ["custom",      "output matches in custom format, i.e. no columns are " +
                  "pre-defined; all columns have to be specified by the user"],
  ["cigar",       "display cigar string representing alignment " +
                  "(no distinction between match and mismatch)"],
  ["cigarX",      "display cigar string representing alignment " +
                  "(distinction between match (=) and mismatch (X))"],
  ["trace",       "display trace, i.e. compact representation of alignment " +
                  "(possibly followed by =<delta> to specify the delta-value, "+
                  "default is #{TRACE_DELTA})"],
  ["dtrace",      "display trace with difference, like trace, but instead " +
                  "of a absolute value x, report #{TRACE_DELTA} - x, which " +
                  "leads to smaller numbers and thus a more compact " +
                  "representation"],
  ["s.len",       "display length of match on subject sequence"],
  ["s.seqnum",    "display sequence number of subject sequence"],
  ["subject id",  "display sequence description of subject sequence"],
  ["s.start",     "display start position of match on subject sequence"],
  ["s.end",       "display end position of match on subject sequence"],
  ["strand",      "display strand of match using symbols F (forward) and " +
                  "P (reverse complement)"],
  ["q.len",       "display length of match on query sequence"],
  ["q.seqnum",    "display sequence number of query sequence"],
  ["query id",    "display sequence description of query sequence"],
  ["q.start",     "display start position of match on query sequence"],
  ["q.end",       "display end position of match on query sequence"],
  ["alignment length",    "display length of alignment"],
  ["mismatches",  "display number of mismatches in alignment"],
  ["indels",      "display number of indels in alignment"],
  ["gap opens",   "display number of indels in alignment"],
  ["score",       "display score of match"],
  ["editdist",    "display unit edit distance"],
  ["identity",    "display percent identity of match"],
  ["seed.len",    "display length seed of the match"],
  ["seed.s",      "display start position of seed in subject"],
  ["seed.q",      "display start position of seed in query"],
  ["s.seqlen",    "display length of subject sequence in which match occurs"],
  ["q.seqlen",    "display length of query sequence in which match occurs"],
  ["evalue",      "display evalue"],
  ["bit score",   "display bit score"]
]

if display_options.length > 64
  STDERR.puts "#{$0}: maximum number of display options is 64"
  exit 1
end

display_options.each do |value|
  if value.length != 2
    STDERR.puts "#{$0}: #{value} is incorrect"
    exit 1
  end
end

kws, max_display_flag_length = keywords(display_options)

outfilename = "src/match/se-display.inc"
begin
  fpout = File.new(outfilename,"w")
rescue => err
  STDERR.puts "#{$0}: cannot create file #{outfilename}"
  exit 1
end

fpout.puts "/* This file was generated by #{$0}, do NOT edit. */"
fpout.puts <<'EOF'
static GtSEdisplayStruct gt_display_arguments_table[] =
{
EOF

kws_sorted = kws.sort {|a,b| a[0] <=> b[0]}
flag2index = Array.new(kws_sorted.length)
kws_sorted.each_with_index do |value,idx|
  flag2index[value[1]] = idx
end

fpout.puts "/* incolumn is true iff the helptext of the argument begins"
fpout.puts "   with the keyword \"display\" */"

fpout.puts kws_sorted.
           map {|s,idx,incolumn| "  {\"#{s}\", #{s.format_enum_value}, #{incolumn}}"}.join(",\n")

fpout.puts <<EOF
};

static unsigned int gt_display_flag2index[] = {
EOF

fpout.puts "   " +  flag2index.join(",\n   ")

fpout.puts <<EOF
};

const char *gt_querymatch_display_help(void)
{
  return "specify what information about the matches to display\\n\"
EOF

longest = 0
display_options.each do |arg,helpline|
  if longest < arg.length
    longest = arg.length
  end
end

display_options.each do |arg,helpline|
  fpout.puts " " * 9 + "\"#{arg}:#{indent(longest,arg)}" +
       "#{format(longest,helpline)}\\n\""
end
fpout.puts <<'EOF'
;
}
EOF

fpout.puts "#define GT_SE_POSSIBLE_DISPLAY_ARGS \"" +
            display_options.map{|arg,helpline| arg + "\""}.
            join("\\\n        \", ")

display_options.each do |arg,helpline|
  if arg == "alignment"
    next
  end
  fpout.puts <<EOF

bool gt_querymatch_#{arg.dot2us}_display(const GtSeedExtendDisplayFlag
                                        *display_flag)
{
  return gt_querymatch_display_on(display_flag,#{arg.format_enum_value});
}
EOF
end

fpout.close_write

outfilename = "src/match/se-display-fwd.inc"
begin
  fpout = File.new(outfilename,"w")
rescue => err
  STDERR.puts "#{$0}: cannot create file #{outfilename}"
  exit 1
end

fpout.puts "/* This file was generated by #{$0}, do NOT edit. */"

fpout.puts "#define GT_DISPLAY_LARGEST_FLAG #{kws.length-1}"
fpout.puts "#define GT_MAX_DISPLAY_FLAG_LENGTH #{max_display_flag_length}"
fpout.puts "#define GT_SEED_EXTEND_DEFAULT_ALIGNMENT_WIDTH #{ALIGNMENT_WIDTH}"
fpout.puts "#define GT_SEED_EXTEND_DEFAULT_TRACE_DELTA #{TRACE_DELTA}"

fpout.puts "typedef enum\n{"

fpout.puts kws.map{|s,idx,incolumn| "  #{s.format_enum_value} /* #{idx} */"}.join(",\n")

fpout.puts "} GtSeedExtendDisplay_enum;"

display_options.each do |arg,helpline|
  if arg != "alignment"
  fpout.puts <<EOF
bool gt_querymatch_#{arg.dot2us}_display(const GtSeedExtendDisplayFlag *);
EOF
  end
end
