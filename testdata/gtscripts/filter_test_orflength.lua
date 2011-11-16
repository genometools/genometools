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

  while not(curnode == nil) do

    if (curnode:get_type() == target) then
      rng = curnode:get_range()
      length = rng:get_end() - rng:get_start() + 1
      if not((length % 3) == 0) then
        return true
      end
    end
    curnode = gfi:next()
  end

  return false
end
