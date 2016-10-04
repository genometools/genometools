#!/usr/bin/env ruby

require 'ostruct'
require 'optparse'

def usage(opts,msg)
  STDERR.puts "#{$0}: #{msg}\n#{opts.to_s}"
  exit 1
end

def parseargs(argv)
  options = OpenStruct.new
  options.typename = :ulong
  ulongpairset = false
  uint64keypairset = false
  flbaset = false
  opts = OptionParser.new
  opts.on("--ulongpair",
          "generate code sorting key/value pairs of ulongs") do |x|
    ulongpairset = true
    options.typename = :ulongkeyvaluepair
  end
  opts.on("--uint64keypair",
          "generate code sorting pairs of uint64_t keys") do |x|
    uint64keypairset = true
    options.typename = :uint64keypair
  end
  opts.on("--flba",
          "generate code sorting byte arrays of fixed length") do |x|
    flbaset = true
    options.typename = :flba
  end
  rest = opts.parse(argv)
  if not rest.empty?
    usage(opts,"superfluous arguments")
  end
  if ulongpairset and uint64keypairset
    STDERR.puts "#{$0}: cannot combine option --ulongpair and --uint64keypair"
    exit 1
  end
  if ulongpairset and flbaset
    STDERR.puts "#{$0}: cannot combine option --ulongpair and --flba"
    exit 1
  end
  if uint64keypairset and flbaset
    STDERR.puts "#{$0}: cannot combine option --uint64keypair and --flba"
    exit 1
  end
  return options
end

def makekeyname(options)
  if options.typename == :ulong
    return "ulong"
  elsif options.typename == :ulongkeyvaluepair
    return "ulongpair"
  elsif options.typename == :uint64keypair
    return "uint64keypair"
  else
    return "flba"
  end
end

def makekeyptr(options)
  return makekeyname(options) + "ptr"
end

def makebasetype(options)
  if options.typename == :ulong
    return "GtUword"
  elsif options.typename == :ulongkeyvaluepair
    return "GtUwordPair"
  elsif options.typename == :uint64keypair
    return "Gtuint64keyPair"
  else
    return "uint8_t"
  end
end

def makevaluetype(options)
  if options.typename == :flba
    return "const " + makebasetype(options) + "*"
  else
    return makebasetype(options)
  end
end

class String
  def del_us
    self.gsub!(/^_/,"")
  end
  def del_us_all
    self.gsub!(/_/,"")
  end
end

def declare_tmpvar(options,var)
  if options.typename == :flba
    return "/* no decl. */"
  else
    return makebasetype(options) + " #{var};"
  end
end

def increment(options,var)
  if options.typename == :flba
    return "#{var} += rbuf->unitsize"
  else
    return "#{var}++"
  end
end

def tmpvarexpr(options,var)
  if options.typename == :flba
    return "rbuf->#{var}_ptr"
  else
    return "#{var}"
  end
end

def offset(options,expr)
  if options.typename == :flba
    return "#{expr} * rbuf->unitsize"
  else
    return expr
  end
end

def copy_ptr_ptr(options,destptr,sourceptr)
  if options.typename == :flba
    return "memcpy(#{destptr},#{sourceptr},rbuf->unitsize)"
  else
    return "*#{destptr} = *#{sourceptr}"
  end
end

def copy_index_expr(options,arr,idx,expr)
  if options.typename == :flba
    if expr.match(/^_/)
      expr = "rbuf->#{expr.del_us}_ptr"
    end
    return "memcpy(#{arr} + (#{idx}) * rbuf->unitsize,\n#{expr},rbuf->unitsize)"
  else
    expr.del_us
    if "#{arr}[#{idx}]".length > 60
      return "#{arr}[#{idx}]=\n#{expr}"
    else
      return "#{arr}[#{idx}] = #{expr}"
    end
  end
end

def copy_var_index(options,var,arr,idx)
  if options.typename == :flba
    if var.match(/^_/)
      var = "rbuf->#{var.del_us}_ptr"
    end
    return "memcpy(#{var},#{arr} + (#{idx}) * rbuf->unitsize,\nrbuf->unitsize)"
  else
    var.del_us
    return "#{var} = #{arr}[#{idx}]"
  end
end

def copy_var_expr(options,var,expr)
  if options.typename == :flba
    if var.match(/^_/)
      var = "rbuf->#{var.del_us}_ptr"
    end
    if expr.match(/^_/)
      expr = "rbuf->#{expr.del_us}_ptr"
    end
    return "memcpy(#{var},#{expr},\nrbuf->unitsize)"
  else
    var.del_us
    expr.del_us
    return "#{var} = #{expr}"
  end
end

def copy_var_ptr(options,var,ptr)
  if options.typename == :flba
    return copy_var_expr(options,var,ptr)
  else
    var.del_us
    return "#{var} = *#{ptr}"
  end
end

def copy_ptr_expr(options,ptr,expr)
  if options.typename == :flba
    return copy_var_expr(options,ptr,expr)
  else
    expr.del_us
    return "*#{ptr} = #{expr}"
  end
end

def shiftnext(options)
  if options.typename == :flba
    return "shift+1"
  else
    return "shift - CHAR_BIT"
  end
end

def recursion_continue(options)
  if options.typename == :flba
    return "currentstackelem.shift < rbuf->unitsize-1"
  else
    return "currentstackelem.shift > 0"
  end
end

def derefptr(options,ptr,comp="a")
  if m = ptr.match(/^&(\w+)/)
    var = m[1]
    var.del_us
  else
    var = ptr
  end
  if options.typename == :ulong
    if var == ptr
      return "*#{var}"
    else
      return "#{var}"
    end
  elsif options.typename == :ulongkeyvaluepair
    if var == ptr
      return "#{var}->#{comp}"
    else
      return "#{var}.#{comp}"
    end
  elsif options.typename == :uint64keypair
    if var == ptr
      return "#{var}->uint64_#{comp}"
    else
      return "#{var}.uint64_#{comp}"
    end
  else
    return ""
  end
end

def radixkey(options,var)
  if options.typename == :ulong or options.typename == :ulongkeyvaluepair
    return "GT_RADIX_KEY(UINT8_MAX,rightshift,#{derefptr(options,var)})"
  elsif options.typename == :uint64keypair
    return "(rightshift > (sizeof (GtUword) - 1) * CHAR_BIT) ?\n" +
           "GT_RADIX_KEY(UINT8_MAX,rightshift - sizeof (GtUword) * CHAR_BIT,\n" +
           "#{derefptr(options,var)}) :\n" +
           "GT_RADIX_KEY(UINT8_MAX,rightshift,#{derefptr(options,var,"b")})"
  else
    if m = var.match(/^&_(\w+)/)
      var = "rbuf->#{m[1]}_ptr"
    end
    return "#{var}[rightshift]"
  end
end

def compare_smaller(options,ptr1,ptr2)
  if options.typename == :ulong or options.typename == :ulongkeyvaluepair
    return "#{derefptr(options,ptr1)} < #{derefptr(options,ptr2)}"
  elsif options.typename == :uint64keypair
    ptr1.del_us_all
    ptr2.del_us_all
    return "gt_radixsort_uint64keypair_smaller(#{ptr1},#{ptr2})"
  else
    if m = ptr1.match(/^&_(\w+)/)
      ptr1 = "rbuf->#{m[1]}_ptr"
    end
    return "memcmp(#{ptr1},#{ptr2},rbuf->unitsize) < 0"
  end
end

options = parseargs(ARGV)

if options.typename == :flba
print <<END_OF_FILE
static #{makebasetype(options)} *gt_radixsort_#{makekeyname(options)}_bin_get(
                                            const GtRadixbuffer *rbuf,
                                            GtUword binnum)
{
  return rbuf->values.#{makekeyptr(options)} +
                 ((binnum << rbuf->log_bufsize) +
                  (GtUword) rbuf->nextidx[binnum]) * rbuf->unitsize;
}
END_OF_FILE
else
print <<END_OF_FILE
static #{makebasetype(options)} gt_radixsort_#{makekeyname(options)}_bin_get(
                                            const GtRadixbuffer *rbuf,
                                            GtUword binnum)
{
  return rbuf->values.#{makekeyptr(options)}[
                 (binnum << rbuf->log_bufsize) +
                 (GtUword) rbuf->nextidx[binnum]];
}
END_OF_FILE
end

print <<END_OF_FILE
static inline void gt_radixsort_#{makekeyname(options)}_bin_update(
                                    #{makebasetype(options)} *source,
                                    GtRadixbuffer *rbuf,
                                    GtUword binnum,
                                    #{makevaluetype(options)} value)
{
  GtUword binoffset = binnum << rbuf->log_bufsize;

  #{copy_index_expr(options,"rbuf->values.#{makekeyptr(options)}\n",
                    "binoffset + (GtUword) rbuf->nextidx[binnum]",
                    "value")};
  if ((GtUword) rbuf->nextidx[binnum] < rbuf->buf_size - 1)
  {
    rbuf->nextidx[binnum]++;
  } else
  {
    GtUword j;
    #{makebasetype(options)} *wsourceptr, *rsourceptr, *rend, *valptr;

    wsourceptr = source +
                 #{offset(options,"(rbuf->endofbin[binnum] - "+
                                  "(rbuf->buf_size - 1))\n")};
    rsourceptr = wsourceptr + #{offset(options,"rbuf->buf_size")};
    rend = source + #{offset(options,"rbuf->startofbin[binnum+1]")};
    valptr = rbuf->values.#{makekeyptr(options)} +
             #{offset(options,"binoffset")};
    for (j=0; j<rbuf->buf_size; j++)
    {
      #{copy_ptr_ptr(options,"wsourceptr","valptr")};
      #{increment(options,"wsourceptr")};
      if (rsourceptr < rend)
      {
        #{copy_ptr_ptr(options,"valptr","rsourceptr")};
        #{increment(options,"rsourceptr")};
      }
      #{increment(options,"valptr")};
    }
    rbuf->nextidx[binnum] = 0;
  }
  rbuf->endofbin[binnum]++;
}

static void gt_radixsort_#{makekeyname(options)}_cached_shuffle(GtRadixbuffer *rbuf,
                                              #{makebasetype(options)} *source,
                                              GtCountbasetype len,
                                              size_t rightshift)
{
  GtUword binoffset, binnum, bufoffset,
                nextbin, firstnonemptybin = UINT8_MAX+1;
  GtCountbasetype *count, previouscount, currentidx;
  #{makebasetype(options)} *sourceptr,
                           *sourceend = source + #{offset(options,"len")};

  rbuf->countcached++;
  count = rbuf->startofbin; /* use same memory for count and startofbin */
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    count[binnum] = 0;
    rbuf->nextidx[binnum] = 0;
  }
  for (sourceptr = source; sourceptr < sourceend; #{increment(options,"sourceptr")})
  {
    count[#{radixkey(options,"sourceptr")}]++;
  }
  for (bufoffset = 0, binoffset = 0, binnum = 0; binnum <= UINT8_MAX;
       bufoffset += rbuf->buf_size, binoffset += count[binnum], binnum++)
  {
    const GtUword elems2copy = MIN(rbuf->buf_size,(GtUword) count[binnum]);

    if (elems2copy > 0)
    {
      if (firstnonemptybin == UINT8_MAX+1)
      {
        firstnonemptybin = binnum;
      }
      memcpy(rbuf->values.
             #{makekeyptr(options)} + #{offset(options,"bufoffset")},
             source + #{offset(options,"binoffset")},
             #{offset(options,"(sizeof *source * elems2copy)")});
    }
  }
  previouscount = count[0];
  rbuf->startofbin[0] = rbuf->endofbin[0] = 0;
  nextbin = 0;
  for (binnum = 1UL; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype temp = rbuf->startofbin[binnum-1] + previouscount;
    previouscount = count[binnum];
    rbuf->startofbin[binnum] = rbuf->endofbin[binnum] = temp;
  }
  /* to simplify compution of bin end */
  rbuf->startofbin[UINT8_MAX+1] = len;
  for (currentidx = 0, binnum = firstnonemptybin;
       currentidx < len; binnum = nextbin - 1)
  {
    #{declare_tmpvar(options,"tmpvalue")}
    #{copy_var_expr(options,"_tmpvalue","gt_radixsort_#{makekeyname(options)}_bin_get(rbuf,binnum)")};
    while (true)
    {
      binnum = #{radixkey(options,"&_tmpvalue")};
      if (currentidx != rbuf->endofbin[binnum])
      {
        #{declare_tmpvar(options,"tmpswap")}
        #{copy_var_expr(options,"_tmpswap","_tmpvalue")};
        #{copy_var_expr(options,"_tmpvalue","gt_radixsort_#{makekeyname(options)}_bin_get(rbuf,binnum)")};
        gt_radixsort_#{makekeyname(options)}_bin_update
                             (source,rbuf,binnum,
                              #{tmpvarexpr(options,"tmpswap")});
      } else
      {
        break;
      }
    }
    gt_radixsort_#{makekeyname(options)}_bin_update(source,rbuf,binnum,
                                           #{tmpvarexpr(options,"tmpvalue")});
    currentidx++;
    /* skip over empty bins */
    while (nextbin <= UINT8_MAX && currentidx >= rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    /* skip over full bins */
    while (nextbin <= UINT8_MAX &&
           rbuf->endofbin[nextbin-1] == rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    if (currentidx < rbuf->endofbin[nextbin-1])
    {
      currentidx = rbuf->endofbin[nextbin-1];
    }
  }
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    GtUword bufleft = (GtUword) rbuf->nextidx[binnum];

    if (bufleft > 0)
    {
      #{makebasetype(options)} *sourceptr, *valptr;

      valptr = rbuf->values.#{makekeyptr(options)} +
               #{offset(options,"(binnum << rbuf->log_bufsize)")};
      sourceptr = source +
                  #{offset(options,"(rbuf->startofbin[binnum+1] - bufleft)")};
      memcpy(sourceptr,valptr,#{offset(options,"(sizeof *sourceptr * bufleft)")});
    }
  }
}

static void gt_radixsort_#{makekeyname(options)}_uncached_shuffle(
                       GtRadixbuffer *rbuf,
                       #{makebasetype(options)} *source,
                       GtCountbasetype len,
                       size_t rightshift)
{
  GtUword binnum, nextbin;
  GtCountbasetype currentidx, previouscount, *count;
  #{makebasetype(options)} *sourceptr,
                           *sourceend = source + #{offset(options,"len")};

  rbuf->countuncached++;
  count = rbuf->startofbin; /* use same memory for count and startofbin */
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    count[binnum] = 0;
    rbuf->nextidx[binnum] = 0;
  }
  for (sourceptr = source; sourceptr < sourceend; #{increment(options,"sourceptr")})
  {
    count[#{radixkey(options,"sourceptr")}]++;
  }
  previouscount = count[0];
  rbuf->startofbin[0] = rbuf->endofbin[0] = 0;
  nextbin = 0;
  for (binnum = 1UL; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype temp = rbuf->startofbin[binnum-1] + previouscount;
    previouscount = count[binnum];
    rbuf->startofbin[binnum] = rbuf->endofbin[binnum] = temp;
  }
  /* to simplify compution of bin end */
  rbuf->startofbin[UINT8_MAX+1] = len;
  for (currentidx = 0; currentidx < len; /* Nothing */)
  {
    GtCountbasetype *binptr;
    #{declare_tmpvar(options,"tmpvalue")}
    #{copy_var_index(options,"_tmpvalue","source","currentidx")};

    while (true)
    {
      binptr = rbuf->endofbin +
               (#{radixkey(options,"&_tmpvalue")});
      binnum = *binptr;
      if (currentidx != binnum)
      {
        #{declare_tmpvar(options,"tmpswap")}
        #{copy_var_expr(options,"_tmpswap","_tmpvalue")};
        #{copy_var_index(options,"_tmpvalue","source","binnum")};
        #{copy_index_expr(options,"source","binnum","_tmpswap")};
        (*binptr)++;
      } else
      {
        break;
      }
    }
    #{copy_index_expr(options,"source","binnum","_tmpvalue")};
    currentidx++;
    (*binptr)++;
    /* skip over empty bins */
    while (nextbin <= UINT8_MAX && currentidx >= rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    /* skip over full bins */
    while (nextbin <= UINT8_MAX &&
           rbuf->endofbin[nextbin-1] == rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    if (currentidx < rbuf->endofbin[nextbin-1])
    {
      currentidx = rbuf->endofbin[nextbin-1];
    }
  }
}

static void gt_radixsort_#{makekeyname(options)}_shuffle(GtRadixbuffer *rbuf,
                                       #{makebasetype(options)} *source,
                                       GtCountbasetype len,
                                       size_t rightshift)
{
  gt_assert(rbuf != NULL);
  if ((GtUword) len > rbuf->cachesize)
  {
    gt_radixsort_#{makekeyname(options)}_cached_shuffle(rbuf,source,len,rightshift);
  } else
  {
    gt_radixsort_#{makekeyname(options)}_uncached_shuffle(rbuf,source,len,
                                                      rightshift);
  }
}

static void
gt_radixsort_#{makekeyname(options)}_inplace_insertionsort(
                                  GT_UNUSED GtRadixbuffer *rbuf,
                                  #{makebasetype(options)} *arr,
                                  GtCountbasetype a_size)
{
  #{makebasetype(options)} *optr,
                           *end = arr + #{offset(options,"a_size")};

  for (optr = arr + #{offset(options,"1")}; optr < end;
       #{increment(options,"optr")})
  {
    #{makebasetype(options)} *oprevious = optr - #{offset(options,"1")};

    if (#{compare_smaller(options,"optr","oprevious")})
    {
      #{makebasetype(options)} *iptr;
      #{declare_tmpvar(options,"tmpvalue")}
      #{copy_var_ptr(options,"_tmpvalue","optr")};

      #{copy_ptr_ptr(options,"optr","oprevious")};
      for (iptr = oprevious; iptr > arr; iptr -= #{offset(options,"1")})
      {
        #{makebasetype(options)} *iprevious = iptr - #{offset(options,"1")};
        if (!(#{compare_smaller(options,"&_tmpvalue","iprevious")}))
        {
          break;
        }
        #{copy_ptr_ptr(options,"iptr","iprevious")};
      }
      #{copy_ptr_expr(options,"iptr","_tmpvalue")};
    }
  }
}

static void gt_radixsort_#{makekeyname(options)}_process_bin(
                                     GtStackGtRadixsort_stackelem *stack,
                                     GtRadixbuffer *rbuf,
                                     #{makebasetype(options)} *source,
                                     size_t shift)
{
  GtUword binnum;

  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype width = rbuf->endofbin[binnum] - rbuf->startofbin[binnum];

    if (width >= (GtCountbasetype) 2)
    {
      #{makebasetype(options)} *ptr
       = source + #{offset(options,"rbuf->startofbin[binnum]")};

      if (width == (GtCountbasetype) 2)
      {
        #{makebasetype(options)} *nextptr = ptr + #{offset(options,1)};
        if (#{compare_smaller(options,"nextptr","ptr")})
        {
          #{declare_tmpvar(options,"tmpswap")}
          #{copy_var_ptr(options,"_tmpswap","ptr")};
          #{copy_ptr_ptr(options,"ptr","nextptr")};
          #{copy_ptr_expr(options,"nextptr","_tmpswap")};
        }
      } else
      {
        if (width <= (GtCountbasetype) 32)
        {
          rbuf->countinsertionsort++;
          gt_radixsort_#{makekeyname(options)}_inplace_insertionsort(rbuf,ptr,width);
        } else
        {
          GtRadixsort_stackelem tmpstackelem;

          tmpstackelem.left.#{makekeyptr(options)} = ptr;
          tmpstackelem.len = width;
          tmpstackelem.shift = #{shiftnext(options)};
          GT_STACK_PUSH(stack,tmpstackelem);
        }
      }
    }
  }
}

static void gt_radixsort_#{makekeyname(options)}_sub_inplace(GtRadixbuffer *rbuf,
                                           GtStackGtRadixsort_stackelem *stack)
{
  GtRadixsort_stackelem currentstackelem;

  while (!GT_STACK_ISEMPTY(stack))
  {
    currentstackelem = GT_STACK_POP(stack);
    gt_radixsort_#{makekeyname(options)}_shuffle(rbuf,
                         currentstackelem.left.#{makekeyptr(options)},
                         currentstackelem.len,
                         currentstackelem.shift);
    if (#{recursion_continue(options)})
    {
      (void) gt_radixsort_#{makekeyname(options)}_process_bin(stack,rbuf,
                                   currentstackelem.left.#{makekeyptr(options)},
                                   currentstackelem.shift);
    }
  }
}
END_OF_FILE
