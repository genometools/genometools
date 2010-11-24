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