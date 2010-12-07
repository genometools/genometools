f = assert(loadfile(arg[1]))
f()

for k,v in pairs(units) do
  if next(v) then
    for innerK, innerV in ipairs(v) do
      io.stdout:write(innerV, ' ')
    end
  end
end
