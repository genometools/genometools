class MultiDimensionalArray
  def initialize(*dimensions)
    @dimensions = Array.new(dimensions.length)
    @factors = Array.new(dimensions.length)
    product = 1
    i = dimensions.length - 1
    while i>=0
      @dimensions[i] = dimensions[i]
      @factors[i] = product
      product *= @dimensions[i]
      i -= 1
    end
    @data = Array.new(product)
  end
  def getOffset(indices)
    if indices.length != @dimensions.length
      STDERR.puts "#{$0}: getOffset indices.length = #{indices.length} " +
                  "!= #{@dimensions.length} = @dimensions.length"
      exit 1
    end
    offset = 0
    0.upto(@dimensions.length-1) do |i|
      if indices[i] < 0 or indices[i] >= @dimensions[i]
        STDERR.puts "#{$0}: illegal value indices[#{i}]=#{indices[i]}"
        exit 1
      end
      offset += @factors[i] * indices[i]
    end
    return offset
  end

  def [](*indices)
    return @data[self.getOffset(indices)]
  end

  def []=(*indicesAndValue)
    value = indicesAndValue.pop
    @data[self.getOffset(indicesAndValue)] = value
  end

  def enumaccesstubles()
    numofwheels = @dimensions.length
    wheelspace = Array.new(numofwheels)
    idx = numofwheels - 1
    0.upto(numofwheels-1) do |i|
      wheelspace[i] = 0
    end
    loop do
      wheelspace[idx] += 1
      if wheelspace[idx] == @dimensions[idx]
        wheelspace[idx] = 0
        if idx == 0
          break
        end
        idx -= 1
      else
        idx = numofwheels-1
        yield wheelspace
      end
    end
  end
end
