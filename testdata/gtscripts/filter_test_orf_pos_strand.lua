name        = "Name here"
author      = "Sascha Kastens"
version     = "1.0"
email       = "sascha.kastens@studium.uni-hamburg.de"
short_descr = "Short description here."
description = "Description here"

function filter(gn)
  target = "reading_frame"
  gfi = gt.feature_node_iterator_new(gn)

  curnode = gfi:next()
  count = 0

  while not(curnode == nil) do

    if (curnode:get_type() == target) then
      if (curnode:get_strand() == '+') then
        count = count + 1
      end
    end
    curnode = gfi:next()
  end

  if (count >= 2) then
    return false
  end

  return true
end
