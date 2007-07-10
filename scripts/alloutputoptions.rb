#!/usr/bin/env ruby

def turnwheels(outlist)
  numofalphabets = outlist.length
  wheelspace = Array.new
  alphasizes = Array.new
  0.upto(numofalphabets-1) do |z|
    alphasizes[z] = 2
    wheelspace[z] = 0
  end
  z = numofalphabets-1
  while true
    output = false
    0.upto(numofalphabets-1) do |i|
      if wheelspace[i] == 1
        output = true
        print " #{outlist[i]}"
      end
    end
    if output
      puts ""
    end
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

turnwheels(["-tis","-suf","-bwt","-lcp"])
