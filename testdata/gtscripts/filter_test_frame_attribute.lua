name        = "ORF filter w/o frame"
author      = "Sascha Kastens"
version     = "1.0"
email       = "sascha.kastens@studium.uni-hamburg.de"
short_descr = "Selects nodes with ORF frame information."
description = "Selects a node if it contains a node of type " ..
              "reading_frame with frame information."

function filter(gn)
  gfi = gt.feature_node_iterator_new(gn)

  target = "reading_frame"

  curnode = gfi:next()
  count = 0

  while not(curnode == nil) do

    if (curnode:get_type() == target) then
      if (curnode:get_attribute("frame") == nil) then
        return true
      end
    end
    curnode = gfi:next()
  end

  return false
end
