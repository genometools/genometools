function usage()
io.stderr:write(string.format("Usage: %s file\n", arg[0]))
  io.stderr:write("Checks a GFF file for line-sortedness.\n")
  os.exit(1)
end

function split(str, sep)
  local fields = {}
  str:gsub("([^"..sep.."]*)"..sep, function(c) table.insert(fields, c) end)
  return fields
end

if #arg == 1 then
  gfffile = arg[1]
else
  usage()
end

cur_seqid = nil
cur_pos = 0
file = assert(io.open(gfffile, "r"))
i = 0
for line in file:lines() do
  i = i + 1
  if not string.match(line, "^#") then
    f = split(line, "\t")
    if #f < 5 then
      io.stderr:write("Not enough fields in line " .. i .. "\n")
      os.exit(1)
    end
    seqid, startpos, endpos = f[1], f[4], f[5]
    if seqid ~= cur_seqid then
      cur_seqid = seqid
      cur_pos = 0
    end
    if tonumber(startpos) < cur_pos then
      io.stderr:write("Error: " .. startpos .. " < "
                       .. cur_pos .. " (line " .. i ..")\n")
      io.stderr:write(line .. "\n")
      os.exit(1)
    end
    cur_pos = tonumber(startpos)
  end
end
file:close()
