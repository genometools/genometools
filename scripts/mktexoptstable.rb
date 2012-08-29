#!/usr/bin/env ruby
# simple transfer -help option list to tex

opts = false
opt = ""
desc = ""
puts <<TEXT
\\begin{table}[hbpt]
  \\centering
  \\caption{<+caption+>}
\\begin{footnotesize}
  \\label{<+label+>}
  \\begin{tabular}{lp{0.6\\textwidth}}
TEXT
ARGF.each_line do |line|
  if line =~ /^-/
    if not opts
      opts = true 
    else
      puts "\\Showoption{#{opt}} & #{desc}\\\\"
    end
    m = line.match /-(\S+)\s+(\S.*)/
    opt = m[1]
    desc = m[2].chomp
  elsif line =~ /^\n/ and opts
    break
  else
    m = line.match /\s+(.*)/
    desc += " " + m[1].chomp
  end
end
puts "\\Showoption{#{opt}} & #{desc}\\\\"
puts <<TEXT
  \\end{tabular}
\\end{footnotesize}
\\end{table}
TEXT
