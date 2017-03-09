def print_sequence(seq,width,fpout = STDOUT)
  idx = 0
  while idx < seq.length
    fpout.puts "#{seq[idx,width]}"
    idx += width
  end
end
