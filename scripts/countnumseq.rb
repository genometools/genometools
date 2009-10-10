def countnumofsequences(inputfile)
  seqcount = 0
  File.open(inputfile).each_line do |line|
    if line.match(/^>/)
      seqcount+=1
    end
  end
  return seqcount
end

