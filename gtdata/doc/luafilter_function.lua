function filter(gn)
  target = "exon"
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
