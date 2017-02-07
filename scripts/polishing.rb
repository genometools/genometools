class String
  def cigar_each(forward = true)
    store = Array.new()
    self.scan(/(\d+)([A-Z])/) do |m|
      multiplier = m[0].to_i
      if not ["M","X","D","I"].member?(m[1])
        raise "cannot parse cigar operation #{m}"
      end
      if forward
        yield multiplier, m[1]
      else
        store.push([multiplier, m[1]])
      end
    end
    if not forward
      store.reverse.each do |multiplier,op|
        yield multiplier, op
      end
    end
  end
end

class Polishing
  def initialize(errorpercentage,history_size)
    if history_size == 0
      @cut_depth = 15;
    else
      @cut_depth = [history_size/2,15].min
    end
    @pol_size = 2 * @cut_depth;
    @match_score = 20.0 * errorpercentage;
    @difference_score = 1000.0 - @match_score;
    puts "match=#{@match_score},difference=-#{@difference_score}"
    @history_size = history_size
  end
  def prefix_positive?(cigarstring)
    return self.end_positive?(cigarstring,true)
  end
  def suffix_positive?(cigarstring)
    return self.end_positive?(cigarstring,false)
  end
  def end_positive?(cigarstring,forward)
    prefix_positive_sum = 0
    hist_rest = @history_size
    cigarstring.cigar_each(forward) do |multiplier,op|
      if multiplier >= hist_rest
        count = hist_rest
        hist_rest = 0
      else
        count = multiplier
        hist_rest -= multiplier
      end
      if op == "M"
        prefix_positive_sum += count * @match_score
      else
        prefix_positive_sum -= count * @difference_score
      end
      if prefix_positive_sum < 0
        return false
      end
      if hist_rest == 0
        break
      end
    end
    return true
  end
end
