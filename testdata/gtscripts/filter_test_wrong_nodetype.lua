name        = "Name here"
author      = "Sascha Kastens"
version     = "1.0"
email       = "mail@skastens.de"
short_descr = "Short description here."
description = "Description here"

function filter(gn)
  target = "xon"
  gfi = gt.feature_node_iterator_new(gn)

  curnode = gfi:next()

  while not(curnode == nil) do

    if (curnode:get_type() == target) then
      return false
    end
    curnode = gfi:next()
  end

  return true
end
