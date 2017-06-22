def turnwheels(asizes)
  numofalphabets = asizes.length
  wheelspace = Array.new
  alphasizes = Array.new
  0.upto(numofalphabets-1) do |z|
    alphasizes[z] = asizes[z]
    wheelspace[z] = 0
  end
  z = numofalphabets-1
  loop do
    yield wheelspace
    stop = false
    while not stop
      wheelspace[z] = wheelspace[z]+1
      if wheelspace[z] == alphasizes[z]
        wheelspace[z] = 0
        if z == 0
          return
        end
        z = z - 1
      else
        z = numofalphabets-1
        stop = true
      end
    end
  end
end
