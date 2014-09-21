math.randomseed(os.time())

function usage()
  io.stderr:write(string.format("Usage: %s indexname minlen maxlen n_substr\n", arg[0]))
  io.stderr:write("Extract <nof_substr> random substrings from a GtEncseq.\n")
  os.exit(1)
end

if #arg == 4 then
  idxname = arg[1]
  minlen = tonumber(arg[2])
  maxlen = tonumber(arg[3])
  nsubstr = tonumber(arg[4])

  el = gt.encseq_loader_new()
  es = el:load(idxname)
  i = 0
  while i < nsubstr do
    len = math.random(minlen, maxlen)
    seqno = math.random(es:num_of_sequences())-1
    eslen = es:seqlength(seqno)
    if eslen > len then
      start = es:seqstartpos(seqno) + math.random(eslen-len)
      stop = start + len - 1
      print(">"..start.."-"..stop.." (length "..len..")")
      print(es:extract_decoded(start, stop))
      i = i + 1
    end
  end
else
  usage()
end
