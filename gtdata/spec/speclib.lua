matchers = {
  should_be = function(value, expected)
    if value ~= expected then
      return false, "expecting "..tostring(expected)..", not ".. tostring(value)
    end
    return true
  end;

  should_be_smaller_than = function(value, expected)
    if value >= expected then
      return false, tostring(value).." is larger than ".. tostring(expected)
    end
    return true
  end;

  should_be_larger_than = function(value, expected)
    if value <= expected then
      return false, tostring(value).." is smaller than ".. tostring(expected)
    end
    return true
  end;

  should_not_be = function(value, expected)
    if value == expected then
      return false, "should not be "..tostring(value)
    end
    return true
  end;

  should_have_key = function(value, expected)
    if value[expected] == nil then
      return false, tostring(value).." does not have key ".. tostring(expected)
    end
    return true
  end;

  should_not_have_key = function(value, expected)
    if value[expected] ~= nil then
      return false, tostring(value).." has key ".. tostring(expected)
    end
    return true
  end;

  should_error = function(f)
    if pcall(f) then
      return false, "expecting an error but received none"
    end
    return true
  end;

  should_match = function(value, pattern)
    if not string.find(value, pattern) then
      return false, value .. " doesn't match pattern "..pattern
    end
    return true
  end;

  should_not_match = function(value, pattern)
    if string.find(value, pattern) then
      return false, value .. " matches pattern "..pattern
    end
    return true
  end;
}
matchers.should_equal = matchers.should_be
