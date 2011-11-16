name        = "LTR node type filter"
author      = "Sascha Kastens"
version     = "1.0"
email       = "sascha.kastens@studium.uni-hamburg.de"
short_descr = "Select nodes with two long terminal repeats."
description = "Selects a node if it contains a node of type " ..
              "LTR_retrotransposon with two child nodes of type " ..
              "long_terminal_repeat."

function filter(gn)
  target = "LTR_retrotransposon"
  subtarget = "long_terminal_repeat"
  gfi = gt.feature_node_iterator_new(gn)

  curnode = gfi:next()
  found_target = false;
  count = 0

  while not(curnode == nil) do
    if (curnode:get_type() == target) then
      found_target = true
    elseif ((found_target == true) and (curnode:get_type() == subtarget)) then
      count = count + 1
    end
    curnode = gfi:next()
  end

  if ((found_target == true) and (count == 2)) then
    return false
  elseif found_target == false then
    return false
  end

  return true
end
